#ifndef AFFINE_HULL_HPP_
#define AFFINE_HULL_HPP_

#include <cassert>

#include <algorithm>
#include <iterator>
#include <vector>

#include <Eigen/Core>

#include "point_cloud.hpp"
#include "dynamic_qr.hpp"

namespace FC {

/**
  @brief supports d + 1 points, where d is the dimension
*/
template <typename point_cloud_t>
class affine_hull {
public:
  typedef point_cloud_t                          point_cloud_type; 
  typedef typename point_cloud_type::number_type number_type;
  typedef typename point_cloud_type::size_type   size_type;
private:
  using member_container = std::vector<size_type>;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using eigen_cmap = Eigen::Map<eigen_vector const>;  // TODO set Aligned flag based on point cloud
public:
  typedef typename member_container::const_iterator iterator;
  
  affine_hull(point_cloud_type const& pc)
    : _dyn_qr(pc.dim()), _members(), _pc(pc) {
  }
  
  template <typename Index>
  void add_point(Index const idx) {
    assert(not is_member(idx));
    assert(size() <= _pc.dim());
    // it requires an idx that serves as the origin for adding a column
    if (size() > 0)
      _dyn_qr.append_column(_pc[idx] - _pc[_members.front()]);
    _members.push_back(idx);
  }
  
  template <typename Index>
  void drop_point(Index const idx) {
    assert(is_member(idx));
    auto const it = std::find(_members.begin(), _members.end(), idx);
    if (_members.size() > 1) {
      bool const del_orig = _members.begin() == it;
      _dyn_qr.delete_column(del_orig ? 0 : // gets deleted because the point
                                           // corresponding to this column
                                           // will become the new origin
                            std::distance(_members.begin(), it) - 1);
      // if there's at least one column left after deleting a member
      // we need to perform a rank-one update on it
      if (del_orig && _members.size() > 2) {
        _dyn_qr.rank_one_update(_pc[_members.front()] - _pc[_members[1]],
                                eigen_vector::Ones(_dyn_qr.num_cols()));
      }
    }
    _members.erase(it);
    assert((_members.empty() && _dyn_qr.num_cols() == 0) ||
           _members.size() == _dyn_qr.num_cols() + 1);
  }
  
  /**
    @return the number of points 
  */
  size_type size() const {
    return _members.size();
  }
  
  iterator begin() const {
    return _members.begin();
  }
  
  iterator end() const {
    return _members.end();
  }
  
  point_cloud_type const& pc() const noexcept {
    return _pc;
  }
  
  /**
    @brief project x onto the affine hull
    @param lambda pointer to the vector of projection coefficients
  */
  template <typename Derived1, typename Derived2>
  void project(Eigen::MatrixBase<Derived1> const& x,
               Eigen::MatrixBase<Derived2> const& lambda) const {
    assert(size() > 0);
    if (_members.size() > 1)
      _dyn_qr.solve(x, const_cast<Eigen::MatrixBase<Derived2> &>(lambda).tail(lambda.size() - 1));
    const_cast<Eigen::MatrixBase<Derived2> &>(lambda)[0]
    = 1 - lambda.tail(lambda.size() - 1).sum();
  }

private:
  bool is_member(size_type const idx) const {
    return _members.end() != std::find(_members.begin(), _members.end(), idx);
  }

  dynamic_qr<number_type> _dyn_qr;
  std::vector<size_type> _members;  // TODO write own dynarray class, to save the capacity pointer
  point_cloud_type const& _pc;
};

}  // namespace FC

#endif  // AFFINE_HULL_HPP_

