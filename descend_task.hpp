#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <Eigen/Dense>

namespace FC {

// forward declaration
template <typename _point_cloud_type>
class ascend_task;

template <typename _point_cloud_type>
class descend_task {
  using number_type = typename _point_cloud_type::number_type;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  
  public:
  
  // TODO almost everything
  
  private:
    eigen_vector _location;
    eigen_vector _ray;

  friend class ascend_task<_point_cloud_type>;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

