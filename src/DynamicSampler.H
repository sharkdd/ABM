#ifndef DynamicSampler_H
#define DynamicSampler_H

#include <iostream>
#include <cassert>
#include <cmath>
#include <list>
#include <numeric>
#include <vector>
//#include <tr1/unordered_map>
#include <unordered_map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <RunningSum.H>

// A data structure for maintaining a discrete probability
// distribution specified by positive sampling weights. The data
// structure supports insertion, deletion, and modification of these
// weights. The implementation is based on Hagerup et al. "Maintaining
// Discrete Probability Distributions Optimally", although it
// sacrifices some computational guarantees (constant time variate
// generation) in favor of flexibility to support floating point
// weights.

// Slepoy et al ("A constant-time kinetic Monte Carlo algorithm for
// simulation of large biochemical reaction networks") presents a more
// accessible discussion of a composition-rejection scheme like that
// implemented here.

// Index must be EqualityComparable
// Hash<Index> must be defined
// Weight is a floating point type
template<typename Index, typename Weight> 
class DynamicSampler {
public:
  // Default constructor. Needs to be initialized separately,
  // otherwise insertion will always fail
  DynamicSampler();

  // Construct with specified minimum and maximum weights on
  // inputs. No separate initialization is required. Weights below
  // wmin may actually be specified, but they are sampled less
  // efficiently.  If the number n of values the distribution may take
  // is fixed, then wmin=1/n^2 is a good choice
  DynamicSampler(Weight const wmin, Weight const wmax);

  // Default destructor
  ~DynamicSampler() {}

  // Initialize the dynamic sampler. If the sampler has previously
  // been initialized, this will force recalculation of the entire
  // data structure.
  void initialize(Weight const wmin, Weight const wmax);

  // Check if an index is in the support of the sampling distribution
  bool member(Index const& index) const {
    return m_offset.find(index) != m_offset.end();
  }

  // Insert an element into the sampler. If the element
  // is already present, its weight is updated
  void insert(Index const& index, Weight const weight);

  // Remove an element from the sampler
  void remove(Index const& index);

  // Update the weight associated with an element
  void update(Index const& index, Weight const weight);

  // Sample an element from the distribution
  Index sample(gsl_rng* rng) const;

  // Return the sum of all stored weights
  Weight sum() const;

  // Return the number of keys
  inline size_t size() const {return m_offset.size();}

private:

  struct Record {
    Record() : weight(0.0) {}
    Record(Index const& i, Weight const w) : index(i), weight(w) {}
    Index index;
    Weight weight;
  };

  struct Offset {
    unsigned int group;
    size_t position;
  };

  typedef std::vector<Record> Group;
  typedef std::unordered_map<Index, Offset> OffsetMap;

  // Returns the group that a weight belongs to
  inline unsigned int group(Weight const weight) {
    int g;
    frexp(weight, &g);
    return g - m_gmin;
  }

  // Remove an index from a group given an iterator to its offset
  // record. This does leave the record in the offset map
  void remove(typename OffsetMap::const_iterator const& it);

  // Helper function used when sampling. Returns the group
  // the sample will be draw from
  unsigned int draw_group(gsl_rng* rng) const;

  // Helper function used when sampling. Returns the position of the
  // sampled index within a selected group
  size_t draw_index(gsl_rng* rng, unsigned int const g) const;

  // Helper function used for sampling from a distribution stored as a
  // vector. This is a template function to allow use with vectors of
  // weights or running sums of weights
  template<typename Value>
  static size_t draw_vector(gsl_rng* rng, std::vector<Value> const& weights);

  OffsetMap m_offset;

  std::vector<Group> m_group;

  // stores group weights
  std::vector< RunningSum<Weight> > m_group_weight;

  // low-weight elements are sampled separately
  std::vector<Weight> m_low_weight;

  RunningSum<Weight> m_sum;

  // upper and lower bounds on the distribution weights
  Weight m_wmin;
  Weight m_wmax;

  // bounds on allowable group indices
  int m_gmin;
  int m_gmax;
};

template<typename Index, typename Weight> 
DynamicSampler<Index,Weight>::DynamicSampler()
  : m_sum(0.0), m_wmin(0.0), m_wmax(0.0) {}

template<typename Index, typename Weight> 
DynamicSampler<Index,Weight>::DynamicSampler(Weight const wmin, Weight const wmax) 
  : m_sum(0.0), m_wmin(wmin), m_wmax(wmax) {
  assert(wmin > 0);
  frexp(wmin, &m_gmin);
  if (frexp(wmax, &m_gmax) > 0.0) ++m_gmax;
  m_group.resize(m_gmax - m_gmin + 1);
  m_group_weight.resize(m_group.size(), 0.0);
}

template<typename Index, typename Weight>
void DynamicSampler<Index,Weight>::initialize(Weight const wmin, Weight const wmax) {
  assert(wmin > 0);

  // Temporary storage for existing records in case of reinitialization
  std::list<Record> records;

  // if the the sampler was already initialized, cache the existing
  // records then clear the data structure
  if (m_wmin != 0 || m_wmax != 0) {
    for (unsigned int g(0); g < m_group.size(); ++g) {
      copy(m_group[g].begin(), m_group[g].end(), records.begin());
      m_group[g].clear();
    }
    m_group_weight.clear();
    m_offset.clear();
  }

  m_wmin = wmin;
  m_wmax = wmax;
  m_sum = 0.0;

  frexp(wmin, &m_gmin);
  if (frexp(wmax, &m_gmax) > 0.0) ++m_gmax;
  m_group.resize(m_gmax - m_gmin + 1);
  m_group_weight.resize(m_group.size(), 0.0);

  // reinsert any records from the previous initialization
  typename std::list<Record>::const_iterator ri;
  for (ri = records.begin(); ri != records.end(); ++ri) insert(ri->index, ri->weight);
}

template<typename Index, typename Weight>
void DynamicSampler<Index,Weight>::insert(Index const& index, Weight const weight) {
  assert(weight <= m_wmax);
  typename OffsetMap::iterator k(m_offset.find(index));
  if (k == m_offset.end()) {
    Offset offset;
    if (weight < m_wmin) {
      offset.group = m_group.size() - 1;
      offset.position = m_low_weight.size();
      m_low_weight.push_back(weight);
    } else {
      offset.group = group(weight);
      offset.position = m_group[offset.group].size();
    }
    m_group[offset.group].push_back(Record(index, weight));
    m_group_weight[offset.group] += weight;
    m_offset.insert(std::pair<Index, Offset>(index, offset));
    m_sum += weight;
  } else {
    update(index, weight);
  }
}

template<typename Index, typename Weight>
void DynamicSampler<Index,Weight>::remove(Index const& index) {
  typename OffsetMap::const_iterator k(m_offset.find(index));
  if (k != m_offset.end()) {
    m_sum -= m_group[k->second.group][k->second.position].weight;
    remove(k);
    m_offset.erase(k);
  }
}

template<typename Index, typename Weight>
void DynamicSampler<Index,Weight>::update(Index const& index, Weight const weight) {
  assert(weight <= m_wmax);
  typename OffsetMap::iterator k(m_offset.find(index));
  if (k != m_offset.end()) {
    unsigned int g;
    g = (weight < m_wmin) ? m_gmin + m_group.size() - 1 : group(weight);
    m_sum += weight - m_group[k->second.group][k->second.position].weight;

    if (g != k->second.group) {
      // weight change changes group membership
      remove(k);
      k->second.group = g;
      k->second.position = m_group[g].size();
      m_group[g].push_back(Record(index, weight));
      m_group_weight[g] += weight;
      if (g == m_group.size() - 1) {
	m_low_weight.push_back(weight);
	assert(m_low_weight.size() - 1 == k->second.position);
      }
    } else {
      // weight change leaves element in the same group
      const size_t j(k->second.position);
      m_group_weight[g] += weight - m_group[g][j].weight;
      m_group[g][j].weight = weight;
      if (g == m_group.size() - 1) m_low_weight[j] = weight;
    }

  } else {
    insert(index, weight);
  }
}

template<typename Index, typename Weight>
void DynamicSampler<Index,Weight>::remove(typename OffsetMap::const_iterator const& it) {
  const unsigned int g(it->second.group);
  if (m_group[g].size() > 1) {
    typename OffsetMap::iterator jt;

    const size_t iold(it->second.position);
    const size_t iend(m_group[g].size() - 1);
    const Weight wold(m_group[g][iold].weight);

    jt = m_offset.find(m_group[g][iend].index);
    assert(jt != m_offset.end());
    jt->second.position = iold;

    m_group_weight[g] -= wold;
    m_group[g][iold] = m_group[g][iend];
    m_group[g].resize(iend);
    if (g == m_group.size() - 1) {
      m_low_weight[iold] = m_low_weight[iend];
    }
  } else {
    m_group_weight[g] = 0.0;
    m_group[g].clear();
    if (g == m_group.size() - 1) m_low_weight.clear();
  }
}

template<typename Index, typename Weight>
unsigned int DynamicSampler<Index,Weight>::draw_group(gsl_rng* rng) const {
  const unsigned int g(draw_vector(rng, m_group_weight));
  assert(m_group[g].size() > 0);
  return g;
}

template<typename Index, typename Weight>
size_t DynamicSampler<Index,Weight>::draw_index(gsl_rng* rng, unsigned int const g) const {
  size_t k;
  if (g == m_group.size() - 1) {
    k = draw_vector(rng, m_low_weight);
  } else {
    // sample among high-weight indices using rejection sampling
    const Weight wmax(gsl_pow_int(2.0, g + m_gmin));
    do k = gsl_rng_uniform_int(rng, m_group[g].size());
    while (gsl_rng_uniform(rng) * wmax > m_group[g][k].weight);
  }
  return k;
}

template<typename Index, typename Weight>
template<typename Value>
size_t DynamicSampler<Index,Weight>::draw_vector(gsl_rng* rng, std::vector<Value> const& weights) {
  RunningSum<Weight> sum(0.0);
  for (size_t k(0); k < weights.size(); ++k) sum += weights[k];

  // gsl_rng_uniform gives us a value in [0,1), so 1-gsl_rng_uniform
  // is in (0,1]. This is because we should never sample an object with zero weight,
  // but if U=0 we will select k=0 if weights[0]=k.
  double U(sum * (1.0 - gsl_rng_uniform(rng)));
  RunningSum<Weight> rsum;
  size_t k(0);
  do {
    rsum += weights[k];
    if (U <= rsum)
      return k;
    ++k;
  } while (k < weights.size());
  // Algebraically, we should never get here so long as
  // sum=sum(weights). The assertion is here in case numeric reality
  // does not meet algebraic expectation.
  assert(false);
  return k;
}

template<typename Index, typename Weight>
Index DynamicSampler<Index,Weight>::sample(gsl_rng* rng) const {
  const unsigned int g(draw_group(rng));
  const size_t k(draw_index(rng, g));
  return m_group[g][k].index;
}

template<typename Index, typename Weight>
Weight DynamicSampler<Index,Weight>::sum() const {
  return static_cast<Weight>(m_sum);
}

#endif // DynamicSampler_H
