#ifndef UPDATE_RAY_HPP_
#define UPDATE_RAY_HPP_

namespace FC {
/**
  @param lambda vector of at least ah.size() elements
*/
template <typename PointCloud, typename Float>
void update_ray(affine_hull<PointCloud> const& ah,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> const& x,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> & lambda,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> & driver,
                Eigen::Matrix<Float, Eigen::Dynamic, 1> & ray) {
  ah.project(x, lambda.head(ah.size()));
  driver.setZero();
  auto m_it = ah.begin();
  for (size_type i = 0; i < ah.size(); ++i)
    driver += lambda[i] * ah.pc()[*m_it++];
  ray = x - driver;
}

}  // namespace FC

#endif  // UPDATE_RAY_HPP_

