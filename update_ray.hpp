#ifndef UPDATE_RAY_HPP_
#define UPDATE_RAY_HPP_

#include <Eigen/Core>

#include "affine_hull.hpp"

namespace FC {

enum class RAY_DIR {FROM_DRIVER, TO_DRIVER};

template <RAY_DIR ray_dir>
struct compute_ray;

/**
  @param lambda vector of at least ah.size() elements
*/
template <RAY_DIR ray_dir, typename PointCloud, typename Float>
void update_ray(affine_hull<PointCloud> const& ah,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> const& x,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> & lambda,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> & driver,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> & ray) {
  using size_type = typename affine_hull<PointCloud>::size_type;
  using eigen_matrix = Eigen::MatrixXd;
  using eigen_vector = Eigen::VectorXd;
  
  ah.project(x, lambda.head(ah.size()));
  driver.setZero();
  auto m_it = ah.begin();
  for (size_type i = 0; i < ah.size(); ++i)
    driver += lambda[i] * ah.pc()[*m_it++];
  compute_ray<ray_dir>::exec(x, driver, ray);
}

template <>
struct compute_ray<RAY_DIR::FROM_DRIVER> {
  template <typename Float>
  static void exec(Eigen::Matrix<Float, Eigen::Dynamic, 1> const& x,
                   Eigen::Matrix<Float, Eigen::Dynamic, 1> const& driver,
                   Eigen::Matrix<Float, Eigen::Dynamic, 1>      & ray) {
  ray = x - driver;
  }
};


template <>
struct compute_ray<RAY_DIR::TO_DRIVER> {
  template <typename Float>
  static void exec(Eigen::Matrix<Float, Eigen::Dynamic, 1> const& x,
                   Eigen::Matrix<Float, Eigen::Dynamic, 1> const& driver,
                   Eigen::Matrix<Float, Eigen::Dynamic, 1>      & ray) {
    ray = driver - x;
  }
};

}  // namespace FC

#endif  // UPDATE_RAY_HPP_

