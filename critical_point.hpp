#ifndef CRITICAL_POINT_HPP_
#define CRITICAL_POINT_HPP_

#include <cassert>

#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>

namespace FC {

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
  critical_point(IdxIterator idx_begin, IdxIterator idx_end, number_type dist)
    : _indices(idx_begin, idx_end), _dist(dist) {}
  
  // constructor for regular cps
  template <typename IdxIterator>
  critical_point(IdxIterator idx_begin, IdxIterator idx_end, number_type dist,
                 self_type * succ)
    : critical_point(idx_begin, idx_end, dist),
      _successors(&succ, &(succ) + 1) {
    assert(succ);
  }
  
  // constructor for cp at inf
  critical_point(size_type index)
    : _index(index) {
  }

  // constructors and assignment-operators
  critical_point(critical_point const&) = delete;
  critical_point(critical_point &&) = default;
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
    assert(_successors.cend() ==
           std::find(_successors.cbegin(), _successors.cend(), succ));
    _successors.push_back(succ);
  }
  
  void erase(self_type *);
  void erase(succ_iterator);

  // info
  bool is_max_at_inf() const noexcept {
    return _indices.empty();
  }
  
  number_type dist() const noexcept {
    return (is_max_at_inf() ? std::numeric_limits<number_type>::infinity()
                            : _dist);
  }
  
  size_type index() const noexcept {
    return (is_max_at_inf() ? _index : std::distance(_indices.cbegin(),
                                                     _indices.cend()));
  }

private:
  idx_container const _indices;
  succ_container      _successors;
  union {
    _number_type const _dist;  // for regular cps
    _size_type const _index;   // for cp at inf
  };
};

template <typename number_type, typename size_type>
bool operator==(critical_point<number_type, size_type> const& lhs,
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

}  // namespace FC

#endif  // CRITICAL_POINT_HPP_
