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

// local headers
#include "compute.hpp"

/**
  @brief This examples is supposed to demonstrate the intended use of the
         library.
*/
int main(int, char**) {
  using number_type = double;
  using size_type = std::uint32_t;
  
  size_type const NUM_PTS  =   10;
  size_type const DIM      =    3;
  size_type const NUM_RUNS = 1000;

  using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
  eigen_matrix data;
  std::array<number_type const*, NUM_PTS> points;

  for (size_type i = 0; i < NUM_RUNS; ++i) {
    unsigned int eigen_seed = std::time(nullptr);
    std::cout << "eigen_seed = " << eigen_seed << std::endl;
    std::srand(eigen_seed);
    
    data = eigen_matrix::Random(DIM, NUM_PTS);
    std::cout << "points = \n" << data << std::endl;
    for (size_type i = 0; i < NUM_PTS; ++i)
      points[i] = data.col(i).data();
    
    auto fc = FC::compute_flow_complex<size_type, false>
              (points.cbegin(), points.cend(), DIM, /* nr of threads */ 8);
    // TODO write dedicated validate function
    std::map<size_type, int> hist;
    for (auto const& cp : fc)
      if (not cp.is_max_at_inf()) {
        auto r_pair = hist.emplace(cp.index(), 1);
        if (not r_pair.second)
          ++r_pair.first->second;
      }
    int sum = 0;
    for (auto const& el : hist)
      sum += (0 == (el.first % 2) ? el.second : -el.second);
    if (1 != sum) {
      std::cout << "HIST-SUM = " << sum << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  std::cout << NUM_RUNS << " RUNS SUCCESSFULLY COMPLETED\n";
  std::exit(EXIT_SUCCESS);
}

