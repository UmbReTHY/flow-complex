#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <utility>

#include <Eigen/Core>

#include "critical_point.hpp"
#include "affine_hull.hpp"
#include "update_ray.hpp"

namespace FC {

template <typename point_cloud_t>
class descend_task {
public:
  typedef point_cloud_t                          point_cloud_type;
  typedef typename point_cloud_type::number_type number_type;
  typedef typename point_cloud_type::size_type   size_type;
private:
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
public:
  descend_task(affine_hull<point_cloud_type> && ah,
               eigen_vector && location,
               critical_point<number_type, size_type> * succ)
    : _ah(std::move(ah)), _succ(succ) {
    _location.swap(location);
  }
  
  template <typename DTHandler, typename ATHandler, typename CPHandler>
  void execute(DTHandler & dth, ATHandler & ath, CPHandler & cph) {
    auto const& pc = _ah.pc();
    eigen_vector ray;
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    update_ray(_ah, _location, lambda, driver, ray);
  }
  
private:
  affine_hull<point_cloud_type>            _ah;
  eigen_vector                             _location;
  critical_point<number_type, size_type> * _succ;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

