#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <Eigen/Dense>

namespace FC {

// forward declaration
template <typename _number_type, typename _size_type>
class ascend_task;

template <typename _number_type, typename _size_type>
class descend_task {
  using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;
  public:
  
  // TODO almost everything
  
  private:
    eigen_vector _location;
    eigen_vector _ray;

  friend class ascend_task<_number_type, _size_type>;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

