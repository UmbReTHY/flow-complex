#ifndef UTILITY_HPP_
#define UTILITY_HPP_

#include <random>
#include <vector>

#include <Eigen/Dense>

namespace FC {

/**
  @brief generates a convex combination of a set of points, where the
         coefficients of the convex combination are strictly greater than 0
  @param begin, end iterators to pointers that hold the coeeficients of the
         points
  @param dim the the size of the point arrays
  @param target a pointer to at least dim elements where the result is placed
         
*/
template <typename number_type, typename PointIterator, typename size_type>
void gen_convex_comb(PointIterator begin, PointIterator end,
                     size_type dim, number_type * target) {
  std::random_device rd;
  std::mt19937 gen(rd());
  // it's important for the lower bound not to be 0, otherwise some points
  // might not take part in the convex combination and the seed point would be
  // located on the convex hull, thus getting stuck during subsequent ascends
  // and not reaching a maximum
  std::uniform_real_distribution<number_type> dis(0.5, 1.0);
  number_type sum = 0.0;
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  Eigen::Map<eigen_vector> target_map(target, dim);
  target_map.setZero();
  number_type tmp;
  for (auto it = begin; it != end; ++it) {
    tmp = dis(gen);
    sum += tmp;
    target_map = tmp * Eigen::Map<eigen_vector const>(*it, dim);
  }
  target_map /= sum;
}

}  // namespace FC

#endif  // UTILITY_HPP_

