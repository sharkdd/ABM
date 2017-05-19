#ifndef BinaryIndexedTree_H
#define BinaryIndexedTree_H

#include <cassert>
#include <vector>

// Binary Indexed Tree implementation, as per
// P. Fenwick, A New Data Structure for Cumulative Frequency Tables,
// Software--Practice and Experience, 1994 24(3):327--336
//
// The type definition from the paper above has been implemented along
// with extensions to support insertion and removal of indices.
//
// The template type T must have group structure with T::operator+
// defined. If t and u have type T, then t+u=u+t, t-t=0, t+0=t;
// additionally T::operator- must be defined such that t+(-u)=t-u.
//
// This data structure may be imprecise if floating point data types
// are used. Consider quantizing values in these cases.
template<typename T> class BinaryIndexedTree {
public:
  BinaryIndexedTree(size_t const n = 0) : m_data(n,0) {}
  ~BinaryIndexedTree() {}

  // Access the value at index k
  T get(size_t k) const;

  // Access the partial sum from index=0 to k
  T sum(size_t k) const;

  // Access the sum of all elements in the tree
  T sum() const {return (size() > 0) ? sum(size()-1) : 0;}

  // Increment the value at index k by increment d
  void inc(size_t k, T d);

  // Set the value at index k equal to value v
  void set(size_t k, T v);

  // Insert the value v and return its index
  size_t insert(T v);

  // Remove the value associated with index k. In general this will
  // rearrange the tree, first by moving the value at the last element
  // to position k, then by reducing the size of the tree by one
  // element.
  size_t remove(size_t k);

  // Find the index with the smallest partial sum that is greater than
  // or equal to v (this may behave differently than expected when
  // v=0). The index returned by BinaryIndexedTree::search(v) is the
  // same as
  //
  // std::lower_bound(V.begin(), V.end(), v) - V.begin()
  //
  // if the partial sums were stored sequentially in the vector V
  size_t search(const T v) const;

  // Return the number of elements stored in the tree
  size_t size() const {return m_data.size();}

  // Empty the tree
  void clear() {m_data.clear();}

private:
  std::vector<T> m_data;
};


template<typename T> T BinaryIndexedTree<T>::get(size_t k) const {
  T v(m_data[k]);
  if (k > 0) {
    size_t parent(k & (k - 1));
    --k;
    while (k != parent) {
      v -= m_data[k];
      k &= k - 1;
    }
  }
  return v;
}


template<typename T> T BinaryIndexedTree<T>::sum(size_t k) const {
  T s(m_data[0]);
  while (k > 0) {
    s += m_data[k];
    k &= k - 1;
  }
  return s;
}

template<typename T> void BinaryIndexedTree<T>::inc(size_t k, T d) {
  const ssize_t K(size()); // suppresses warnings when comparing i to size()
  ssize_t i(k); // need signed (two's complement) i for (i & -i) to work
  do {
    m_data[i] += d;
    i += (i & -i);
  } while ((i > 0) && (i < K));
}

template<typename T> void BinaryIndexedTree<T>::set(size_t k, T v) {
  inc(k, v - get(k));
}

template<typename T> size_t BinaryIndexedTree<T>::insert(T v) {
  size_t k(m_data.size());
  m_data.push_back(0);
  m_data[k] = -get(k);
  inc(k, v);
  return k;
}

template<typename T> size_t BinaryIndexedTree<T>::remove(size_t k) {
  assert(k < size());
  size_t n(size() - 1);
  inc(k, get(n) - get(k));
  m_data.resize(n);
  return n;
}

template<typename T>
size_t BinaryIndexedTree<T>::search(const T v) const {
  size_t kmin(0), kmax(size()-1), kmid;
  if (v <= sum(kmin)) return kmin;

  // It is somewhat more effective to use an initial pivot that is a
  // power of two, since fewer locations are accessed to evaluate the
  // partial sum at these indices.

  kmid = 1;
  do kmid <<= 1; while (kmid < size());
  kmid >>= 1;

  // standard binary search
  while (kmax - kmin > 1) {
    if (v > sum(kmid)) kmin = kmid;
    else kmax = kmid;
    kmid = (kmax + kmin) >> 1;
  }
  return kmax;
}

#endif // BinaryIndexedTree_H
