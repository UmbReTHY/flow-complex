#ifndef AFFINE_HULL_HPP_
#define AFFINE_HULL_HPP_

#include <cassert>

#include <algorithm>
#include <iterator>
#include <vector>

#include <Eigen/Dense>

#include "point_cloud.hpp"
#include "dynamic_qr.hpp"

namespace FC {

/**
  @brief supports d + 1 points, where d is the dimension
*/
template <typename _number_type, typename _size_type>
class affine_hull {
  using member_container = std::vector<_size_type>;
  using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
  using eigen_cmap = Eigen::Map<eigen_vector const>;

  public:
    typedef _number_type                                 number_type;
    typedef _size_type                                     size_type;
    typedef typename member_container::const_iterator const_iterator;
    typedef point_cloud<number_type, size_type>        point_cloud_t;
    
    affine_hull(point_cloud_t const& pc)
      : _dyn_qr(pc.dim()), _members(), _pc(pc) {
    }
    
    void add_point(size_type const idx) {
      assert(!is_member(idx));
      assert(size() <= _pc.dim());
      // it requires an idx that serves as the origin for adding a column
      if (size() > 0)
        _dyn_qr.append_column(eigen_cmap(_pc[idx], _pc.dim()) -
                              eigen_cmap(_pc[_members.front()], _pc.dim()));
      _members.push_back(idx);
    }
    
    void drop_point(size_type const idx) {
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
          _dyn_qr.rank_one_update(eigen_cmap(_pc[_members.front()], _pc.dim()) -
                                  eigen_cmap(_pc[_members[1]], _pc.dim()),
                                  eigen_vector::Ones(_dyn_qr.num_cols()));
        }
      }
      _members.erase(it);
      assert((_members.empty() && _dyn_qr.num_cols() == size_type(0)) ||
             _members.size() == _dyn_qr.num_cols() + size_type(1));
    }
    
    /**
      @return the number of points 
    */
    size_type size() const {
      return _members.size();
    }
    
    const_iterator begin() const {
      return _members.begin();
    }
    
    const_iterator end() const {
      return _members.end();
    }
    
    /**
      @brief project x onto the affine hull
      @param lambda pointer to the vector of projection coefficients
    */
    void project(number_type const* x, number_type * lambda) {
      assert(size() > 0);
      if (_members.size() > 1)
        _dyn_qr.solve(x, lambda + 1);
      *lambda = 1 - eigen_cmap(lambda + 1, _dyn_qr.num_cols()).sum();
    }
  
  private:
    bool is_member(_size_type const idx) const {
      return _members.end() != std::find(_members.begin(), _members.end(), idx);
    }

    dynamic_qr<_number_type, _size_type> _dyn_qr;
    std::vector<_size_type> _members;  // TODO write own dynarray class, to save the capacity pointer
    point_cloud_t const& _pc;
};

}  // namespace FC

#endif  // AFFINE_HULL_HPP_

