#ifndef VERTEX_FILTER_HPP_
#define VERTEX_FILTER_HPP_

#include <cassert>

#include "affine_hull.hpp"
#include "point_cloud.hpp"

namespace FC {

template <typename PointCloud>
struct vertex_filter {
  typedef typename PointCloud::size_type size_type;
  typedef typename PointCloud::eigen_map eigen_map;
  vertex_filter(affine_hull<PointCloud> const& ah) : _ah(ah), _curr_idx(0) {
  }
  
  void reset() {
    _curr_idx = 0;
  }

  eigen_map const* operator()(size_type * idx_ptr) {
    assert(idx_ptr);
    while (_ah.end() != std::find(_ah.begin(), _ah.end(), _curr_idx))
      ++_curr_idx;  // skip affine hull members
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
};

// helper function
template <typename PointCloud>
vertex_filter<PointCloud>
make_vertex_filter(affine_hull<PointCloud> const& ah) {
  return vertex_filter<PointCloud>(ah);
}

}  // namespace FC

#endif  // VERTEX_FILTER_HPP_

