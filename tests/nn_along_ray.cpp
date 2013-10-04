
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ctime>

#include <vector>

#include <Eigen/Dense>

#include "nn_along_ray.hpp"
#include "point_cloud.hpp"

using number_type = long double;
using size_type = std::size_t;
using dim_type = std::uint8_t;

dim_type const DIM = 8;
size_type const NUM_PTS = 50;
size_type const NUM_MEMBERS = 3;
static_assert(NUM_MEMBERS < NUM_PTS, "");
// how often a random 
size_type const NUM_TRIALS = 1000;

using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;

using namespace fc;

eigen_matrix points;
eigen_vector loc;
std::vector<size_type> members;

void generate_data() {
  number_type const EMPTY_BALL_RADIUS = number_type(1.0);
  number_type const MAX_RADIUS = number_type(10.0);
  static_assert(MAX_RADIUS > EMPTY_BALL_RADIUS, "");
  eigen_matrix points(DIM, NUM_PTS);
  std::vector<size_type> members;
  // generate members at random on sphere
  size_type i;
  for (i = 0; i < NUM_MEMBERS; ++i) {
    members.push_back(i);
    eigen_vector const tmp = eigen_vector::Random(DIM);
    points.col(i) = tmp * (EMPTY_BALL_RADIUS / tmp.norm());
  }
  
  // generate rest of pts at random outside the same sphere
    for (; i < NUM_PTS; ++i) {
      eigen_vector const tmp = eigen_vector::Random(DIM);
      auto const rnd_num = std::rand() / static_cast<number_type>(RAND_MAX);
      // all other pts should be outside the member sphere
      auto const small_offset = 1e-5;
      auto const radius = (EMPTY_BALL_RADIUS + small_offset) +
                          (MAX_RADIUS - EMPTY_BALL_RADIUS) * rnd_num;
      points.col(i) = tmp * (radius / tmp.norm());
  }
  
  // select the location and shift all points to that origin
  loc = eigen_vector::Random(DIM);
  points += loc.replicate(1, NUM_PTS);
}

int main() {
  std::srand(std::time(nullptr));
  generate_data();
  std::vector<number_type const*> ptr_cont;
  for (size_type i = 0; i < pt_src.cols(); ++i)
    ptr_cont.push_back(pt_src.col(i).data());

  point_cloud<number_type, size_type, dim_type> pc(DIM, ptr_cont.begin(),
                                                        ptr_cont.end());
  
  for (size_type i = 0; i < NUM_TRIALS; ++i) {
    eigen_vector ray = eigen_vector::Random(DIM);
    auto const res_pair =
    nearest_neighbor_along_ray(loc.data(), ray.data(), );
    generate_data();
  }
  
  std::exit(EXIT_SUCCESS);
}
