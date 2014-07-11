#ifndef VERTEX_FILTER_HPP_
#define VERTEX_FILTER_HPP_

#include <cassert>

#include <algorithm>
#include <ostream>

#include "affine_hull.hpp"
#include "logger.hpp"

namespace FC {

/**
  @brief a vertex filter optimized specifically for the conditions of a
         descend_task
*/
template <class PointCloud, class ResultIterator>
struct vertex_filter {
  typedef typename PointCloud::size_type size_type;
  typedef typename PointCloud::eigen_map eigen_map;

  /**
    @brief for descend tasks
  */
  template <class Derived, class MemberIterator>
  vertex_filter(PointCloud const& pc,
                Eigen::MatrixBase<Derived> const& driver,
                MemberIterator member_begin, MemberIterator member_end,
                size_type ignore_idx, ResultIterator result_begin)
  : _pc(pc), _current(result_begin), _end(result_begin) {
    _end = pc.radius_search(driver,
                            (pc[*member_begin] - driver).squaredNorm(),
                            _current);
    _end = std::remove_if(_current, _end, [=](size_type idx) {
      return (member_end != std::find(member_begin, member_end, idx)) or
             idx == ignore_idx;
    });
    Logger() << "VERTEX-FILTER ctor: candidate indices = ";
    for (auto it = _current; it != _end; ++it)
      Logger() << *it << ", ";
    Logger() << std::endl;
  }
  
  /**
    @brief for ascend tasks
  */
  template <class Derived1, class Derived2, class MemberIterator>
  vertex_filter(PointCloud const& pc,
                Eigen::MatrixBase<Derived1> const& ah_member,
                Eigen::MatrixBase<Derived2> const& v,
                MemberIterator member_begin, MemberIterator member_end,
                size_type ignore_idx, ResultIterator result_begin)
  : _pc(pc), _current(result_begin), _end(result_begin) {
    at_init(ah_member, v, member_begin, member_end, ignore_idx);
    Logger() << "VERTEX-FILTER ctor: candidate indices = ";
    for (auto it = _current; it != _end; ++it)
      Logger() << *it << ", ";
    Logger() << std::endl;
  }
  
  /**
    @brief only called by ascend tasks
  */
  template <class Derived>
  void reset(Eigen::MatrixBase<Derived> const& ray,
             affine_hull<PointCloud> const& ah, ResultIterator result_begin) {
    Logger() << "reset vertex-filter\n";
    _current = result_begin; _end = result_begin;
    at_init(_pc[*ah.begin()], ray, ah.begin(), ah.end(), *ah.begin());
  }
  
  eigen_map const* operator()(size_type * idx_ptr) {
    assert(idx_ptr);
    eigen_map const* r;
    if (_current != _end) {
      *idx_ptr = *_current++;
      r = &_pc[*idx_ptr];
    } else {
      r = nullptr;
    }
    return r;
  }
  
private:
  template <class Derived1, class Derived2, class MemberIterator>
  void at_init(Eigen::MatrixBase<Derived1> const& ah_member,
               Eigen::MatrixBase<Derived2> const& v,
               MemberIterator member_begin, MemberIterator member_end,
               size_type ignore_idx) {
    using float_t = typename PointCloud::number_type;
    // precomputed values
    float_t const vdota = v.dot(ah_member);
    for (size_type i = 0; i < _pc.size(); ++i) {
      if ((v.dot(_pc[i]) - vdota) > 0 and
          member_end == std::find(member_begin, member_end, i) and
          i != ignore_idx)
        *(_end++) = i;
    }
  }

  PointCloud const& _pc;
  ResultIterator    _current;
  ResultIterator    _end;
};

// helper function
template <typename PointCloud, class Iterator, class Derived>
vertex_filter<PointCloud, Iterator>
make_at_filter(affine_hull<PointCloud> const& ah,
               Eigen::MatrixBase<Derived> const& ray,
               typename PointCloud::size_type ignore_idx,
               Iterator result_begin) {
  using vf_type = vertex_filter<PointCloud, Iterator>;
  return vf_type(ah.pc(), ah.pc()[*ah.begin()], ray,
                 ah.begin(), ah.end(), ignore_idx, result_begin);
}

template <class PointCloud, class Iterator, class Derived>
vertex_filter<PointCloud, Iterator>
make_dt_filter(affine_hull<PointCloud> const& ah,
               Eigen::MatrixBase<Derived> const& driver,
               typename PointCloud::size_type ignore_idx,
               Iterator result_begin) {
  using vf_type = vertex_filter<PointCloud, Iterator>;
  return vf_type(ah.pc(), driver, ah.begin(), ah.end(), ignore_idx,
                 result_begin);
}

}  // namespace FC

#endif  // VERTEX_FILTER_HPP_

