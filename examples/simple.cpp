// C headers
#include <cstdint>
#include <cstdlib>
#include <ctime>

// C++ headers
#include <array>
#include <iostream>

#define EIGEN_DONT_VECTORIZE

// 3rd-party library headers
#include <Eigen/Core>

// local headers
#include "compute.hpp"

/**
  @brief This example is supposed to demonstrate the intended use of the
         library.
*/
int main(int, char**) {
  // parameters of the point cloud and the computation
  using size_type = std::uint32_t;
  size_type constexpr NUM_PTS  =     200;
  size_type constexpr DIM      =       3;
  size_type constexpr NUM_THREADS =    1;
  // create random point cloud
  using eigen_matrix = Eigen::MatrixXd;
  using float_t = Eigen::MatrixXd::Scalar;
  std::srand(std::time(nullptr));
  eigen_matrix const point_storage = eigen_matrix::Random(DIM, NUM_PTS);
  std::array<float_t const*, NUM_PTS> points;
  for (size_type i = 0; i < NUM_PTS; ++i)
    points[i] = point_storage.col(i).data();
  // compute the flow complex
  auto fc = FC::compute_flow_complex<size_type>
            (points.cbegin(), points.cend(), DIM, NUM_THREADS);

  if (not FC::validate(fc)) {
    std::cerr << "flow complex validation failed\n";
    std::exit(EXIT_FAILURE);
  }
  std::cout << "FC COMPUTATION SUCCESSFULLY COMPLETED\n";
  std::exit(EXIT_SUCCESS);
}

