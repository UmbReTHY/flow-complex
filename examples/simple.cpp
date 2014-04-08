// C headers
#include <cstdint>
#include <cstdio>

// C++ headers
#include <array>

// 3rd-party library headers
#include <Eigen/Core>

// local headers
#include "critical_point.hpp"
#include "compute.hpp"

template <typename number_type, typename size_type>
void print(FC::critical_point<number_type, size_type> const&);

/**
  @brief This examples is supposed to demonstrate the intended use of the
         library.
*/
int main(int, char**) {
  using number_type = double;
  using size_type = std::uint32_t;
  
  size_type const NUM_PTS = 1000;
  size_type const DIM     =    3;

  // create a point cloud - uniform 1000 pt sampling of DIM-d space
  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  eigen_matrix data = eigen_matrix::Random(DIM, NUM_PTS);
  std::array<number_type const*, NUM_PTS> points;
  for (size_type i = 0; i < NUM_PTS; ++i)
    points[i] = data.col(i).data();

  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using cmap = Eigen::Map<eigen_vector const>;

  auto fc = FC::compute_flow_complex<size_type, false>
            (points.cbegin(), points.cend(), DIM, 
             /* nr of threads */ 8, /* numerical tolerance */ 1.e-5);
  
  for (auto const& cp : fc)
    print(cp);
}

template <typename number_type, typename size_type>
void print(FC::critical_point<number_type, size_type> const& cp) {
  using std::printf;
  printf("index of cp = %u\n", cp.index());
  if (!cp.is_max_at_inf()) {
    printf("\tpoint-indices: ");
    for (auto it = cp.idx_cbegin(); it != cp.idx_cend(); ++it)
      printf("%u, ", *it);
    printf("\n");
    printf("\tnumber of successors: %ld\n",
           std::distance(cp.succ_cbegin(), cp.succ_cend()));
    printf("\ttheir distances: ");
    for (auto it = cp.succ_cbegin(); it != cp.succ_cend(); ++it)
      printf("%f, ", (*it)->dist());
    printf("\n");
  }
  printf("\tdist = %f\n", cp.dist());
}
