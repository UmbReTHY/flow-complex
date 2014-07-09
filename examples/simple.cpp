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
#include "logger.hpp"

/**
  @brief This examples is supposed to demonstrate the intended use of the
         library.
*/
int main(int, char**) {
  using number_type = double;
  using size_type = std::uint32_t;
  
  size_type const NUM_PTS  =  500;
  size_type const DIM      =    3;
  size_type const NUM_RUNS =    1;

  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;

  std::array<eigen_vector, NUM_PTS>       point_storage;
  std::array<number_type const*, NUM_PTS>        points;
  // init raw (aligned) memory and point storage
  for (size_type i = 0; i < NUM_PTS; ++i) {
    point_storage[i] = eigen_vector(DIM);
    points[i] = point_storage[i].data();
    if (((unsigned long)points[i] & 15) != 0) {
      std::cout << "misaligned\n";
      std::exit(EXIT_FAILURE);
    }
  }

  for (size_type i = 0; i < NUM_RUNS; ++i) {
    unsigned int eigen_seed = 1568751282;// std::time(nullptr);
    FC::Logger() << "eigen_seed = " << eigen_seed << std::endl;
    std::srand(eigen_seed);
    
    FC::Logger() << "points = \n:";
    for (size_type i = 0; i < NUM_PTS; ++i) {
      point_storage[i].setRandom();
      FC::Logger() << point_storage[i].transpose() << std::endl;
    }
    
    auto fc = FC::compute_flow_complex<size_type, true>
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
    for (auto const& el : hist) {
      std::cout << "index " << el.first << " : " << el.second << std::endl;
      sum += (0 == (el.first % 2) ? el.second : -el.second);
    }
    if (1 != sum) {
      std::cout << "HIST-SUM = " << sum << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  std::cout << NUM_RUNS << " RUNS SUCCESSFULLY COMPLETED\n";
  std::exit(EXIT_SUCCESS);
}

