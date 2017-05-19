#ifndef DynamicSampler_H
#define DynamicSampler_H

#include <stdint.h>
#include <cassert>
#include <cmath>
#include <limits>
#include <list>
#include <numeric>
#include <vector>
#include <tr1/unordered_map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>


// This dynamic sampler implementation uses quantized weights. Weights
// are stored as 32-bit unsigned integers. The largest representable
// weight W depends on the Scale template parameter, with
// W=2^(32-log2(Scale))-1. Keys with weight below 1/Scale will not be
// sampled
//
// This is designed to hold at most 2^32-1 unique keys with weights in
// the range [0,2^16-1). 
template<typename Key, uint16_t Scale = 65535>
class DynamicSampler {
private:
  typedef uint32_t Weight;
  typedef unsigned long long sum_type;

  struct Offset {uint16_t group; uint32_t index;};

  typedef std::tr1::unordered_map<Key, Offset> KeyMap;

  struct Record {
    Record() {}
    Record(Key const& k, Weight const w) : key(k), weight(w) {}
    Key key;
    Weight weight;
  };

  typedef std::vector<Record> Group;

public:

  DynamicSampler();

  ~DynamicSampler() {
    delete [] m_group;
    delete [] m_group_weight;
  }

  inline size_t bucket_count() const {return m_keymap.bucket_count();}
  inline void rehash(size_t const n) {m_keymap.rehash(n);}

  // Check if an index is in the support of the sampling distribution
  inline bool member(Key const& key) const {
    return m_keymap.find(key) != m_keymap.end();
  }

  // Insert an element into the sampler. If the element is already
  // present, its weight is updated
  void insert(Key const& key, double const weight);

  // Remove an element from the sampler
  void remove(Key const& key);

  // Update the weight associated with a key
  inline void update(Key const& key, double const weight) {insert(key, weight);}

  // Sample a key from the distribution
  Key sample(gsl_rng* rng) const;

  // Return the sum of all stored weights
  inline double sum() const {return m_sum / static_cast<double>(Scale);}

  // Return the number of keys
  inline size_t size() const {return m_keymap.size();}

private:

  static uint16_t weight2group(Weight weight) {
    uint16_t value(1);
    while (weight) {++value; weight >>= 1;}
    return value;
  }

  static inline Weight maximum(uint16_t const group) {
    return 1 << (group - 1);
  }

  static inline Weight quantize(double const weight) {
    return static_cast<Weight>(floor(Scale * weight));
  }

  // Helper function that handles modifying the weight of a key
  // already present in the sampler
  void update(typename KeyMap::iterator const& iterator, Weight const weight);

  // helper function to remove a key from a group. This does not
  // remove the keymap record for that key. The result of this
  // operation is undefined if the group or index is out of bounds
  inline void remove_from_group(uint16_t const group, uint32_t const index) {
    const size_t n(m_group[group].size());
    m_group[group][index] = m_group[group][n - 1];
    m_keymap[m_group[group][index].key].index = index;
    m_group[group].resize(n - 1);
  }

  // Helper functions used for sampling the group and an index within
  // a group, respectively
  uint16_t draw_group(gsl_rng* rng) const;
  uint32_t draw_index(gsl_rng* rng, uint16_t const group) const;

  double m_max_weight;
  uint16_t m_num_groups;
  KeyMap m_keymap;
  Group *m_group;
  sum_type *m_group_weight;
  sum_type m_sum;
};

template<typename Key, uint16_t Scale>
DynamicSampler<Key, Scale>::DynamicSampler()
  : m_max_weight(std::numeric_limits<uint32_t>::max() / static_cast<double>(Scale)),
    m_sum(0) {
  m_num_groups = weight2group(quantize(m_max_weight)) + 1;
  m_group = new Group[m_num_groups];
  m_group_weight = new sum_type[m_num_groups];
  for (uint16_t group(0); group < m_num_groups; ++group) {
    m_group_weight[group] = 0;
  }
}


template<typename Key, uint16_t Scale>
void DynamicSampler<Key, Scale>::insert(Key const& key, double const weight) {
  assert((weight >= 0.0) && (weight < m_max_weight));

  const Weight w(quantize(weight));
  typename KeyMap::iterator ki(m_keymap.find(key));
  if (ki == m_keymap.end()) {
    // insert the new element
    Offset offset;
    Record record(key, w);

    offset.group = weight2group(w);
    offset.index = m_group[offset.group].size();
    m_keymap.insert(std::pair<Key,Offset>(key, offset));

    m_group[offset.group].push_back(record);
    m_group_weight[offset.group] += w;
    m_sum += w;

  } else {
    update(ki, w);
  }
}


template<typename Key, uint16_t Scale>
void DynamicSampler<Key, Scale>::remove(Key const& key) {
  typename KeyMap::iterator ki(m_keymap.find(key));
  if (ki != m_keymap.end()) {
    const uint16_t group(ki->second.group);
    const uint32_t index(ki->second.index);
    const Weight weight(m_group[group][index].weight);

    remove_from_group(group, index);
    m_keymap.erase(ki);

    m_group_weight[group] -= weight;
    m_sum -= weight;
  }
}


template<typename Key, uint16_t Scale>
void DynamicSampler<Key, Scale>::update(typename KeyMap::iterator const& iterator, Weight const w_new) {
  const uint16_t g_old(iterator->second.group), g_new(weight2group(w_new));
  const uint32_t i_old(iterator->second.index);
  const Weight w_old(m_group[g_old][i_old].weight);

  m_sum -= w_old;
  m_sum += w_new;
  m_group_weight[g_old] -= w_old;
  if (g_new == g_old) { // the key remains in its current group
    m_group[g_new][i_old].weight = w_new;   // update the key's weight
    m_group_weight[g_new] += w_new;
  } else { // the key changes group
    const size_t n_new(m_group[g_new].size());

    // overwrite the updated key's position in its old group using the
    // last key in the group, update bookkeeping for the latter key,
    // then shrink the group
    remove_from_group(g_old, i_old);

    // place the key in its new group
    m_group[g_new].push_back(Record(iterator->first, w_new));

    // update the offset stored in the key map
    iterator->second.group = g_new;
    iterator->second.index = n_new;

    // update the new group's weight
    m_group_weight[g_new] += w_new;
  }
}


template<typename Key, uint16_t Scale>
Key DynamicSampler<Key, Scale>::sample(gsl_rng* rng) const {
  const uint16_t group(draw_group(rng));
  const uint32_t index(draw_index(rng, group));
  return m_group[group][index].key;
}


template<typename Key, uint16_t Scale>
uint16_t DynamicSampler<Key, Scale>::draw_group(gsl_rng* rng) const {
  // Get a random variable U uniformly distributed in (0, m_sum]
  const double U(m_sum * (1.0 - gsl_rng_uniform(rng)));
  uint16_t group;
  sum_type running_sum(0);
  for (group = 0; group < m_num_groups; ++group) {
    running_sum += m_group_weight[group];
    if (U <= running_sum) return group;
  }
  assert(false); // Should never get here
  return group;
}


template<typename Key, uint16_t Scale>
uint32_t DynamicSampler<Key, Scale>::draw_index(gsl_rng* rng, uint16_t const group) const {
  const double max_weight(maximum(group));
  uint32_t index;
  do index = gsl_rng_uniform_int(rng, m_group[group].size());
  while (gsl_rng_uniform(rng) * max_weight > m_group[group][index].weight);
  return index;
}

#endif // DynamicSampler_H
