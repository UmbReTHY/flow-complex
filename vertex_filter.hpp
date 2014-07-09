#ifndef VERTEX_FILTER_HPP_
#define VERTEX_FILTER_HPP_

#include <cassert>

#include <algorithm>
#include <ostream>

#include "affine_hull.hpp"
#include "logger.hpp"

namespace FC {

template <typename PointCloud, class Iterator>
struct vertex_filter {
  typedef typename PointCloud::size_type size_type;
  typedef typename PointCloud::eigen_map eigen_map;
  vertex_filter(affine_hull<PointCloud> const& ah,
                Iterator drop_begin, Iterator drop_end)
    : _ah(ah), _curr_idx(0), _drop_begin(drop_begin), _drop_end(drop_end) {
    Logger() << "VERTEX-FILTER ctor: dropped-indices = ";
    for (auto it = _drop_begin; it != _drop_end; ++it)
      Logger() << *it << ", ";
    Logger() << std::endl;
  }
  
  void reset(Iterator drop_end) {
    _curr_idx = 0;
    _drop_end = drop_end;
    Logger() << "VERTEX-FILTER reset(): dropped-indices = ";
    for (auto it = _drop_begin; it != _drop_end; ++it)
      Logger() << *it << ", ";
    Logger() << std::endl;
  }

  eigen_map const* operator()(size_type * idx_ptr) {
    assert(idx_ptr);
    while (_ah.end() != std::find(_ah.begin(), _ah.end(), _curr_idx) or
           _drop_end != std::find(_drop_begin, _drop_end, _curr_idx))
      ++_curr_idx;  // skip affine hull members, and lately dropped indices
    eigen_map const* r;
    if (_curr_idx < _ah.pc().size()) { // valid index
      *idx_ptr = _curr_idx;
      r = & _ah.pc()[_curr_idx++];  // this ++ is important!!!
    } else {
      r = nullptr;
    }
    return r;
  }
  
private:
  affine_hull<PointCloud> const& _ah;
  size_type                      _curr_idx;
  Iterator                       _drop_begin;
  Iterator                       _drop_end;
};

/**
  @brief a vertex filter optimized specifically for the conditions of a
         descend_task
*/
template <class PointCloud, class ResultIterator>
struct dt_vertex_filter {
  typedef typename PointCloud::size_type size_type;
  typedef typename PointCloud::eigen_map eigen_map;

  template <class Derived, class MemberIterator>
  dt_vertex_filter(PointCloud const& pc,
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
  PointCloud const& _pc;
  ResultIterator    _current;
  ResultIterator    _end;
};

// helper function
template <typename PointCloud, class Iterator>
vertex_filter<PointCloud, Iterator>
make_vertex_filter(affine_hull<PointCloud> const& ah,
                   Iterator drop_begin, Iterator drop_end) {
  return vertex_filter<PointCloud, Iterator>(ah, drop_begin, drop_end);
}

template <class PointCloud, class Iterator, class Derived>
dt_vertex_filter<PointCloud, Iterator>
make_vertex_filter(affine_hull<PointCloud> const& ah,
                   Eigen::MatrixBase<Derived> const& driver,
                   typename PointCloud::size_type ignore_idx,
                   Iterator result_begin) {
  using vf_type = dt_vertex_filter<PointCloud, Iterator>;
  return vf_type(ah.pc(), driver, ah.begin(), ah.end(), ignore_idx,
                 result_begin);
}

}  // namespace FC

#endif  // VERTEX_FILTER_HPP_

