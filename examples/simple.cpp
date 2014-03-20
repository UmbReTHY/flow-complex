// C headers
#include <cstdint>

// C++ headers
#include <iostream>
#include <ostream>
#include <vector>
#include <memory>

// 3rd-party library headers
#include <Eigen/Dense>

// local headers
#include "critical_point.hpp"
#include "compute.hpp"
#include "flow_complex.hpp"

template <typename number_type, typename size_type>
std::ostream & operator<<(std::ostream &,
                          FC::critical_point<number_type, size_type> const&);

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
  std::vector<number_type const*> points(NUM_PTS);
  for (size_type i = 0; i < data.cols(); ++i)
    points[i] = data.col(i).data();

  auto fc = FC::compute_flow_complex(points.data(), points.data() + points.size(),
                                     DIM, 
                                     /* nr of threads */ 8,
                                     /* numerical tolerance */ 1.e-5);
  
  for (auto const& cp : fc)
    std::cout << cp << std::endl;
}

template <typename number_type, typename size_type>
std::ostream & operator<<(std::ostream & out_stream,
                          FC::critical_point<number_type, size_type> const& cp) {
  out_stream << "index of cp = " << cp.index() << std::endl;
  if (!cp.is_max_at_inf()) {
    out_stream << "\tpoint-indices: ";
    for (auto it = cp.idx_cbegin(); it != cp.idx_cend(); ++it)
      out_stream << *it << ", ";
    out_stream << std::endl;
    
    out_stream << "\tnumber of successors: " << std::distance(cp.succ_cbegin(), cp.succ_cend()) << std::endl;
    out_stream << "\ttheir distances: ";
    for (auto it = cp.succ_cbegin(); it != cp.succ_cend(); ++it)
      out_stream << (*it)->dist() << ", ";
    out_stream << std::endl;
  }
  out_stream << "\tdist = " << cp.dist();
  
  return out_stream;
}
