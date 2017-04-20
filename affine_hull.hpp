#ifndef AFFINE_HULL_HPP_
#define AFFINE_HULL_HPP_

#include <cassert>

#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "dynamic_qr.hpp"
#include "utility.hpp"

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
public:
  typedef typename member_container::const_iterator const_iterator;
  
  affine_hull(point_cloud_type const& pc)
    : _dyn_qr(pc.dim()), _members(), _pc(pc) {
    DLOG(INFO) << "**AH-CTOR " << this << std::endl;
  }
  
  affine_hull(affine_hull const& orig)
    : _dyn_qr(orig._dyn_qr), _members(orig._members), _pc(orig._pc) {
    DLOG(INFO) << "**AH-COPY-CTOR " << this << std::endl;
  }
  
  affine_hull(affine_hull && tmp)
    : _dyn_qr(std::move(tmp._dyn_qr)), _members(std::move(tmp._members)),
      _pc(tmp._pc) {
    DLOG(INFO) << "**AH-MOVE-CTOR " << this << std::endl;
  }

  ~affine_hull() {
    DLOG(INFO) << "**AH-DESTRUCT " << this << std::endl;
  }

  affine_hull & operator=(affine_hull &&) = delete;
  affine_hull & operator=(affine_hull const&) = delete;
  
  void append_point(size_type const idx) {
    DCHECK(_members.end() == std::find(_members.begin(), _members.end(), idx))
    << "added a point to the affine hull that was already contained: " << idx;
    DCHECK(size() <= _pc.dim()) << "added more than d+1 points to affine hull";
    // it requires an idx that serves as the origin for adding a column
    if (size() > 0)
      _dyn_qr.append_column(_pc[idx] - _pc[_members.front()]);
    _members.push_back(idx);
  }
  
  /**
    @brief invalidates all iterators at and after it
  */
  void drop_point(const_iterator it) {
    DLOG(INFO) << "dropping point " << *it;
    DCHECK(is_member(it));
    if (_members.size() > 1) {
      bool const del_orig = _members.begin() == it;
      _dyn_qr.delete_column(del_orig ? 0 : // gets deleted because the point
                                           // corresponding to this column
                                           // will become the new origin
                            std::distance(_members.cbegin(), it) - 1);
      // if there's at least one column left after deleting a member
      // we need to perform a rank-one update on it
      if (del_orig && _members.size() > 2) {
        _dyn_qr.rank_one_update(_pc[_members.front()] - _pc[_members[1]],
                                eigen_vector::Ones(_dyn_qr.num_cols()));
      }
    }
    // instead of using std::vector::erase, we use this workaround to both
    // preserve the order without the element being deleted and support
    // standard library implementations before C++11
    DCHECK(!_members.empty());
    const std::size_t del_pos = std::distance(_members.cbegin(), it);
    auto mutable_it = std::next(_members.begin(), del_pos);
    std::rotate(mutable_it, std::next(mutable_it), _members.end());
    _members.resize(_members.size() - 1);
    DCHECK((_members.empty() && _dyn_qr.num_cols() == 0) ||
           (_members.size() == _dyn_qr.num_cols() + 1));
  }
  
  /**
    @return the number of points 
  */
  size_type size() const {
    return convertSafelyTo<size_type>(_members.size());
  }
  
  const_iterator begin() const {
    return _members.begin();
  }
  
  const_iterator end() const {
    return _members.end();
  }
  
  const point_cloud_type & pc() const noexcept {
    return _pc;
  }
  
  /**
    @brief project x onto the affine hull
    @param lambda pointer to the vector of projection coefficients
  */
  template <typename Derived1, typename Derived2>
  void project(Eigen::MatrixBase<Derived1> const& x,
               Eigen::MatrixBase<Derived2> const& lambda_const) const {
    assert(size() > 0);
    assert(lambda_const.size() == size());
    using lambda_t = Eigen::MatrixBase<Derived2>;
    auto & lambda = const_cast<lambda_t &>(lambda_const);
    if (_members.size() > 1)
      _dyn_qr.solve(x - _pc[*begin()], lambda.tail(lambda_const.size() - 1));
    lambda[0] = 1 - lambda_const.tail(lambda_const.size() - 1).sum();
  }

private:
  bool is_member(const_iterator it) const {
    return (_members.begin() <= it) && (it < _members.end());
  }

  dynamic_qr<number_type> _dyn_qr;
  member_container _members;
  point_cloud_type const& _pc;
};

template <class PointCloud>
std::ostream & operator<<(std::ostream & os, affine_hull<PointCloud> const& ah) {
  os << "HULL MEMBERS: ";
  for (auto it = ah.begin(); it != ah.end(); ++it)
    os << *it << ", ";
  return os;
}

}  // namespace FC

#endif  // AFFINE_HULL_HPP_

