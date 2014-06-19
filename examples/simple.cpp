// C headers
#include <cstdint>
#include <cstdlib>
#include <ctime>

// C++ headers
#include <array>
#include <iostream>
#include <map>

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
  
  size_type const NUM_PTS = 15;
  size_type const DIM     =    3;

  // create a point cloud - uniform 1000 pt sampling of DIM-d space
  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  
  unsigned int eigen_seed = 1403214468;//std::time(nullptr);
  std::srand(eigen_seed);
  
  eigen_matrix data = eigen_matrix::Random(DIM, NUM_PTS);
  std::array<number_type const*, NUM_PTS> points;
  for (size_type i = 0; i < NUM_PTS; ++i)
    points[i] = data.col(i).data();

  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using cmap = Eigen::Map<eigen_vector const>;
  
  int sum;
//  do {
    std::cout << "eigen_seed = " << eigen_seed << std::endl;
    auto fc = FC::compute_flow_complex<size_type, false>
              (points.cbegin(), points.cend(), DIM, /* nr of threads */ 8);
    
    std::map<size_type, int> hist;
    for (auto const& cp : fc)
      if (not cp.is_max_at_inf()) {
        auto r_pair = hist.emplace(cp.index(), 1);
        if (not r_pair.second)
          ++r_pair.first->second;
      }
    sum = 0;
    for (auto const& el : hist)
      sum += (0 == (el.first % 2) ? el.second : -el.second);
    std::cout << "HIST-SUM = " << sum << std::endl;
//  } while (sum == 1);
//    print(cp);
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
      if ((*it)->is_max_at_inf())
        cout << "inf, ";
      else
        cout << (*it)->sq_dist() << ", ";
    }
    cout << '\n';
    cout << "\tsquared dist = " << cp.sq_dist();
  } else {
    cout << "\tdist = inf";
  }
  cout << std::endl;
}
