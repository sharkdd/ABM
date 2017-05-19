#ifndef DYNAMIC_QUEUE_H
#define DYNAMIC_QUEUE_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <tr1/unordered_map>


// A queue of (key, value) pairs prioritized by value. Keys are
// unique; pushing a (key, value) pair that has the same key as a
// previous pair will update the value of the existing pair.
//
// K is an arbitrary key type
// V are the values used to prioritize entries in the queue
// C is a comparator with operator() defined, such that if C(v,u)=true
// then u has higher priority than v
template <typename K, typename V, typename C = std::less<V> >
class DynamicQueue {
private:
  // The DynamicQueue must support efficient bidirectional mapping,
  // since updating a key's priority requires finding it's position in
  // the queue, and we also need to find the highest priority key
  // quickly. To keep time and memory requirements consistent, this is
  // managed by (1) using a hash table (KeyToIndex) to map keys to
  // queue indices and (2) storing priorities and pointers to hash
  // table values in the queue (in seperate vectors, both
  // heap-ordered). The heap stores pointers (as opposed to hash table
  // iterators) because unordered_map rehashing invalidates iterators
  // but not pointers.
  typedef size_t Index;
  typedef typename std::tr1::unordered_map<K,Index> KeyToIndex;
  typedef typename KeyToIndex::value_type* P;
  typedef typename std::vector<P> IndexToKey;
  typedef typename std::vector<V> IndexToVal;
public:
  typedef typename IndexToKey::size_type size_type;
  typedef typename std::pair<K,V> value_type;

  DynamicQueue() {}
  ~DynamicQueue() {}

  size_type size() const {return m_key.size();}
  void clear() {m_ind.clear(); m_key.clear(); m_val.clear();}
  bool empty() const {return size() == 0;}

  bool member(K const& key) const {return m_ind.find(key) != m_ind.end();}

  // returns the number of elements erased (0 or 1)
  size_type erase(K const& key);

  const value_type top() const {return value_type(m_key[0]->first, m_val[0]);}
  void pop();
  void push(K const& key, V const& val) {insert(key,val);}

  void dbg_print(std::ostream &out) {
    for (size_type k(0); k < size(); ++k) {
      out << m_key[k]->second << ' '
	  << m_key[k]->first << ' '
	  << m_val[k] << '\n';
    }
  }

private:
  void insert(K const& key, V const& val);
  void heap_rise(Index index);
  void heap_fall(Index index);

  KeyToIndex m_ind; // maps each key to a unique index
  IndexToKey m_key; // maps each index to a key, heap-ordered by C::operator()
  IndexToVal m_val; // maps each index to its value, heap-ordered by C::operator()

  C m_comp; // comparator
};

template <typename K, typename V, typename C>
void DynamicQueue<K,V,C>::insert(K const& key, V const& val) {
  typename KeyToIndex::iterator it(m_ind.find(key));
  if (it == m_ind.end()) {
    it = m_ind.insert(it, std::pair<K,Index>(key, size()));
    m_key.push_back(&(*it));
    m_val.push_back(val);
  } else {
    m_val[it->second] = val;
    heap_fall(it->second);
  }
  return heap_rise(it->second);
}

template <typename K, typename V, typename C>
typename DynamicQueue<K,V,C>::size_type DynamicQueue<K,V,C>::erase(K const& key) {
  typename KeyToIndex::iterator it(m_ind.find(key));

  if (it == m_ind.end()) return 0;

  // replace with the rightmost element 
  const size_type k(it->second), n(size()-1);
  if (k < n) {
    m_key[k] = m_key[n];
    m_val[k] = m_val[n];
    m_key.resize(n);
    m_val.resize(n);
    m_ind.erase(it);

    // keep the hash table synchronized with the heap
    m_key[k]->second = k;

    // move the replaced element up/down the heap
    heap_rise(it->second);
    heap_fall(it->second);
  } else {
    m_key.resize(n);
    m_val.resize(n);
    m_ind.erase(it);
  }

  return 1;
}

template <typename K, typename V, typename C>
void DynamicQueue<K,V,C>::pop() {
  typename KeyToIndex::iterator it;
  it = m_ind.find(m_key[0]->first);

  const size_type n(size()-1);
  m_ind.erase(it);
  if (n > 0) {
    m_key[0] = m_key[n];
    m_val[0] = m_val[n];
    m_key.resize(n);
    m_val.resize(n);

    m_key[0]->second = 0;
    heap_fall(0);
  } else {
    m_key.resize(n);
    m_val.resize(n);
  }
}

template <typename K, typename V, typename C>
void DynamicQueue<K,V,C>::heap_rise(Index n) {
  P p(m_key[n]);
  for (size_type k(n); k > 0; n = k) {
    k = (n - 1) >> 1;
    if (m_comp(m_val[k], m_val[n])) {
      m_key[k]->second = n;
      std::swap(m_key[k], m_key[n]);
      std::swap(m_val[k], m_val[n]);
    } else break;
  }
  p->second = n;
}

template <typename K, typename V, typename C>
void DynamicQueue<K,V,C>::heap_fall(Index n) {
  P p(m_key[n]);
  Index k, m;
  V v(m_val[n]), vk, vm;
  for (;; n = k) {
    k = (n << 1) + 1;
    m = (n << 1) + 2;
    vk = (k < size() ? m_val[k] : v);
    vm = (m < size() ? m_val[m] : v);
    if (m_comp(v, vk)) {
      if (m_comp(vk, vm)) k = m;
    } else if (m_comp(v, vm)) {
      k = m;
    } else break;
    m_key[k]->second = n;
    std::swap(m_key[k], m_key[n]);
    std::swap(m_val[k], m_val[n]);
  }
  p->second = n;
}

#endif // DYNAMIC_QUEUE_H
