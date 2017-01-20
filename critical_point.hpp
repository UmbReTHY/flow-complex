#ifndef CRITICAL_POINT_HPP_
#define CRITICAL_POINT_HPP_

#include <cassert>

#include <algorithm>
#include <iterator>
#include <vector>
#include <ostream>
#include <istream>

#include <glog/logging.h>
#include <tbb/mutex.h>

#include "utility.hpp"
#include "makros.h"

namespace FC {

// forward declaration
template <typename nt, typename st> class flow_complex;

/**
  @brief 
*/
template <typename _number_type, typename _size_type>
class critical_point {
public:
  typedef _number_type                           number_type;
  typedef _size_type                               size_type;
  typedef critical_point<_number_type, _size_type> self_type;
private:
  using succ_container = std::vector<self_type *>;
  using idx_container = std::vector<_size_type>;
public:    
  typedef typename idx_container::const_iterator   idx_iterator;
  typedef typename succ_container::const_iterator succ_iterator;

  // constructor for finite maxima
  template <typename IdxIterator>
  critical_point(IdxIterator idx_begin, IdxIterator idx_end,
                 number_type sq_dist)
    : _indices(idx_begin, idx_end), _sq_dist(sq_dist) {
    DLOG(INFO) << "**CP-CTOR: " << this << std::endl;
    std::sort(_indices.begin(), _indices.end());
  }
  
  // constructor for regular cps
  template <typename IdxIterator>
  critical_point(IdxIterator idx_begin, IdxIterator idx_end,
                 number_type sq_dist, self_type * succ)
    : critical_point(idx_begin, idx_end, std::move(sq_dist)) {
    DLOG(INFO) << "**CP-CTOR: " << this << std::endl;
    assert(succ);
    _successors.push_back(succ);
  }
  
  // constructor for cp at inf
  critical_point(size_type index)
    : _index(index) {
    DLOG(INFO) << "**CP-INF-CTOR: " << this << std::endl;
  }

  // constructors and assignment-operators
  critical_point(critical_point && tmp) : _indices(std::move(tmp._indices)),
    _successors(std::move(tmp._successors)) {
    DLOG(INFO) << "**CP-MOVE-CTOR: " << this << std::endl;
    if (is_max_at_inf())
      _index = std::move(tmp._index);
    else
      _sq_dist = std::move(tmp._sq_dist);
  }
  
  critical_point(critical_point const& orig)
    : _indices(orig._indices), _successors(orig._successors) {
    if (orig.is_max_at_inf())
      _index = orig._index;
    else
      _sq_dist = orig._sq_dist;
  }

  critical_point & operator=(critical_point const&) = delete;
  critical_point & operator=(critical_point &&) = delete;

  // idx-iterators
  idx_iterator idx_begin() const noexcept {
    return _indices.cbegin();
  }
  
  idx_iterator idx_end() const noexcept {
    return _indices.cend();
  }

  // succ iterators
  succ_iterator succ_begin() const noexcept {
    return _successors.cbegin();
  }
  
  succ_iterator succ_end() const noexcept {
    return _successors.cend();
  }
  
  // modifiers
  void add_successor(self_type * succ) {
    tbb::mutex::scoped_lock lock(_succ_mutex);
    if (succ_end() == std::find(succ_begin(), succ_end(), succ))
      _successors.push_back(succ);
  }
  
  // this method is used by the cleansing tool only, and so far not in a multi-
  // threaded fashion (therefore no mutex overhead).
  // Instead of using std::vector::erase right away, this workaround is supposed
  // to handle both standard library implementations (before and after C++11)
  void erase(self_type const* succ) {
    auto it = std::find(_successors.begin(), _successors.end(), succ);
    if (it != succ_end()) {
      std::swap(*it, _successors.back());
      _successors.resize(_successors.size() - 1);
    }
  }
  void erase(succ_iterator);

  // info
  bool is_max_at_inf() const noexcept {
    return _indices.empty();
  }
  
  /**
    @return if is_max_at_inf() is true, then the result is undefined
  */
  number_type sq_dist() const noexcept {
    return _sq_dist;
  }
  
  size_type index() const noexcept {
    return (is_max_at_inf() ? _index : std::distance(_indices.cbegin(),
                                                     _indices.cend()) - 1);
  }

private:
  idx_container  _indices;
  succ_container _successors;
  tbb::mutex     _succ_mutex;
  union {
    _number_type _sq_dist;  // for regular cps
    _size_type   _index;   // for cp at inf
  };
  
  template <typename nt, typename st>
  friend std::istream & operator>>(std::istream &, flow_complex<nt, st> &);
};

template <typename number_type, typename size_type>
__pure bool operator==(critical_point<number_type, size_type> const& lhs,
                       critical_point<number_type, size_type> const& rhs) {
  assert(std::is_sorted(lhs.idx_begin(), lhs.idx_end()));
  assert(std::is_sorted(rhs.idx_begin(), rhs.idx_end()));
  // we want to allow cps of different indices to be compared
  // (e.g. when stored in the same hash container, == is possibly called on
  // two instances with different index), however this should evaluate to
  // false then
  return ((lhs.is_max_at_inf() && rhs.is_max_at_inf()) ||     // both cps at inf
          ((!lhs.is_max_at_inf() && !rhs.is_max_at_inf()) &&  // both regular
           (lhs.index() == rhs.index()) &&  // avoids comparison below (somet.)
            std::equal(lhs.idx_begin(), lhs.idx_end(), rhs.idx_begin())
           )
         );
}

template <typename number_type, typename size_type>
bool operator!=(critical_point<number_type, size_type> const& lhs,
                critical_point<number_type, size_type> const& rhs) {
  return !(lhs == rhs);
}

struct CPHash {
  template <class CriticalPoint>
  std::size_t operator()(CriticalPoint const& cp) const {
    RangeHash range_hash;
    return range_hash(cp.idx_begin(), cp.idx_end());
  }
};

/**
  @brief prints only the index sequence
*/
template <typename nt, typename st>
std::ostream & operator<<(std::ostream & os, critical_point<nt, st> const& cp) {
  if (cp.is_max_at_inf()) {
    os << "inf ";
  } else {
    for (auto it = cp.idx_begin(); it != cp.idx_end(); ++it)
      os << *it << " ";
  }
  return os;
}

}  // namespace FC

#endif  // CRITICAL_POINT_HPP_
