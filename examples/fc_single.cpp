// C headers
#include <cstdint>
#include <cstdlib>
#include <ctime>

// C++ headers
#include <array>
#include <iostream>
#include <chrono>

#define EIGEN_DONT_VECTORIZE

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
  
  size_type constexpr NUM_PTS  =  120;
  size_type constexpr DIM      =    4;
  size_type constexpr NUM_RUNS =   10;

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

  float time_sum = 0;
  for (size_type i = 0; i < NUM_RUNS; ++i) {
    unsigned int eigen_seed = 1568751282;// std::time(nullptr);
    FC::Logger() << "eigen_seed = " << eigen_seed << std::endl;
    std::srand(eigen_seed);
    
    FC::Logger() << "points = \n:";
    for (size_type i = 0; i < NUM_PTS; ++i) {
      point_storage[i].setRandom();
      FC::Logger() << point_storage[i].transpose() << std::endl;
    }
    
    using sysclock = std::chrono::system_clock;
    using time_pt = sysclock::time_point;
    time_pt start = sysclock::now();
    auto fc = FC::compute_flow_complex<size_type, true>
             (points.cbegin(), points.cend(), DIM, 1);
    time_pt end = sysclock::now();
    using namespace std::chrono;
    time_sum += duration_cast<milliseconds>(end - start).count();

    if (not FC::validate(fc)) {
      std::cerr << "flow complex validation failed\n";
      std::exit(EXIT_FAILURE);
    }
  }
  std::cout << "single thread time : " << (time_sum / NUM_RUNS) << "ms\n";
  
  std::cout << NUM_RUNS << " RUNS SUCCESSFULLY COMPLETED\n";
  std::exit(EXIT_SUCCESS);
}

