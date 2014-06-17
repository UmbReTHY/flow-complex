// C headers
#include <cstdint>

// C++ headers
#include <array>
#include <iostream>

// 3rd-party library headers
#include <Eigen/Core>
//#include <boost/multiprecision/gmp.hpp> 

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
//  namespace mp = boost::multiprecision;
//  using number_type = mp::number<mp::gmp_float<256>, mp::et_off>;
  using number_type = double;
  using size_type = std::uint32_t;
  
  size_type const NUM_PTS = 1000;
  size_type const DIM     =    3;
//  size_type const NUM_PTS = 4;
//  size_type const DIM     =    2;

  // create a point cloud - uniform 1000 pt sampling of DIM-d space
  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  eigen_matrix data = eigen_matrix::Random(DIM, NUM_PTS);
//  data << 1.0, 0.5, 1.0, 6.0,
//          0.25, 1.5, 2.75, 1.5;
  std::array<number_type const*, NUM_PTS> points;
  for (size_type i = 0; i < NUM_PTS; ++i)
    points[i] = data.col(i).data();

  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using cmap = Eigen::Map<eigen_vector const>;

  auto fc = FC::compute_flow_complex<size_type, false>
            (points.cbegin(), points.cend(), DIM, /* nr of threads */ 8);
  
  for (auto const& cp : fc)
    print(cp);
}

template <typename number_type, typename size_type>  // TODO summary (alternating sum)
void print(FC::critical_point<number_type, size_type> const& cp) {
  using std::cout;
  cout << "index of cp = " << cp.index() << '\n';
  if (!cp.is_max_at_inf()) {
    cout << "\tpoint-indices: ";
    for (auto it = cp.idx_begin(); it != cp.idx_end(); ++it)
      cout << *it << ", ";
    cout << '\n';
    cout << "\tnumber of successors: "
         << std::distance(cp.succ_begin(), cp.succ_end()) << '\n';
    cout << "\ttheir squared distances: ";
    for (auto it = cp.succ_begin(); it != cp.succ_end(); ++it) {
      cout << (*it)->sq_dist() << ", ";
    }
    cout << '\n';
    cout << "\tsquared dist = " << cp.sq_dist();
  } else {
    cout << "\tdist = inf";
  }
  cout << std::endl;
}
