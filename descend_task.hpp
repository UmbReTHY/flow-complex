#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <utility>

#include <Eigen/Dense>

#include "affine_hull.hpp"

namespace FC {

// forward declaration
template <typename point_cloud_t>
class ascend_task;

template <typename point_cloud_t>
class descend_task {
  using point_cloud_type = point_cloud_t;
  using number_type = typename point_cloud_type::number_type;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  
public:
  typedef point_cloud_t point_cloud_type;

  descend_task(affine_hull<point_cloud_type> && ah,
               eigen_vector && location,
               eigen_vector && ray) : _ah(std::move(ah)) {
    _location.swap(location);
    _ray.swap(ray);
  }
  
  template <typename DTHandler, typename ATHandler, typename CPHandler>
  void execute(DTHandler & dth, ATHandler & ath, CPHandler & cph) {
  }
  
private:
  affine_hull<point_cloud_type> _ah;
  eigen_vector                  _location;
  eigen_vector                  _ray;
  // TODO remove this: better provide a ascend_task-constructor
  //                   that directly takes the required arguments
  friend class ascend_task<point_cloud_type>;  
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

