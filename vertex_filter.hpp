#ifndef VERTEX_FILTER_HPP_
#define VERTEX_FILTER_HPP_

#include <cassert>

#include "affine_hull.hpp"
#include "point_cloud.hpp"

namespace FC {

template <typename PointCloud, class Iterator>
struct vertex_filter {
  typedef typename PointCloud::size_type size_type;
  typedef typename PointCloud::eigen_map eigen_map;
  // TODO: respect dropped indices
  vertex_filter(affine_hull<PointCloud> const& ah,
                Iterator drop_begin, Iterator drop_end)
    : _ah(ah), _curr_idx(0), _drop_begin(drop_begin), _drop_end(drop_end) {
  }
  
  void reset(Iterator drop_end) {
    _curr_idx = 0;
    _drop_end = drop_end;
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

// helper function
template <typename PointCloud, class Iterator>
vertex_filter<PointCloud, Iterator>
make_vertex_filter(affine_hull<PointCloud> const& ah,
                   Iterator drop_begin, Iterator drop_end) {
  return vertex_filter<PointCloud, Iterator>(ah, drop_begin, drop_end);
}

}  // namespace FC

#endif  // VERTEX_FILTER_HPP_

