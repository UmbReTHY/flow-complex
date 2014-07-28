#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "nn_along_ray.hpp"
#include "point_cloud.hpp"

using number_type = double;
using size_type = std::size_t;
using dim_type = std::uint8_t;
using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;

using res_t = std::pair<number_type, std::vector<size_type>>;

// TODO use CATCH for the test cases

int main() {
  // variables needed through the entire test
  dim_type DIM = 2;
  size_type NUM_PTS;
  eigen_matrix data;
  eigen_vector loc;
  eigen_vector ray;
  std::vector<number_type const*> ptr_cont;
  std::vector<size_type> members;
  res_t res;
  using namespace fc;
  using pc_t = point_cloud<number_type, size_type, dim_type>;

  // 1) no nn found
  NUM_PTS = 7;
  data = eigen_matrix(DIM, NUM_PTS);
  data << 0.5, 1.0, 1.5, 2.5, 3.5, 4.0, 4.5,
          0.5, 1.5, 1.0, 0.5, 1.0, 1.5, 0.5;
  ptr_cont.clear();
  for (size_type i = 0; i < data.cols(); ++i)
    ptr_cont.push_back(data.col(i).data());
  pc_t pc1 (DIM, ptr_cont.begin(), ptr_cont.end());
  members.clear();
  members = {1, 5};
  loc = eigen_vector(DIM);
  loc << 2.5, 2.5;
  ray = eigen_vector(DIM);
  ray << 0.0, 1.0;
  res = nearest_neighbor_along_ray(loc.data(), ray.data(), pc1,
                                   members.begin(), members.end());
  assert(res.second.empty());
  
  // 2) 2 nn found, after 2 had to be discarded
  NUM_PTS = 11;
  data = eigen_matrix(DIM, NUM_PTS);
  data << 0.5, 1.5, 1.5, 2.5, 3.5, 3.5, 4.5, 2.0, 3.0, 1.5, 3.5,
          0.5, 2.0, 1.0, 0.5, 1.0, 2.0, 0.5, 5.5, 5.5, 5.0, 5.0;
  ptr_cont.clear();
  for (size_type i = 0; i < data.cols(); ++i)
    ptr_cont.push_back(data.col(i).data());
  pc_t pc2 (DIM, ptr_cont.begin(), ptr_cont.end());
  res = nearest_neighbor_along_ray(loc.data(), ray.data(), pc2,
                                   members.begin(), members.end());
  auto const& nn = res.second;
  std::vector<size_type> const expected_nn = {9, 10};
  assert(nn.size() == 2 &&
         std::is_permutation(nn.begin(), nn.end(), expected_nn.begin()));
  assert(Eigen::internal::isApprox(res.first, 1.0));
  
  std::exit(EXIT_SUCCESS);
}

