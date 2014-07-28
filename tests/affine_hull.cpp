#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ctime>

#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>

#include <Eigen/Dense>

#include "affine_hull.hpp"

// TODO convert this to CATCH test framework

using number_type = double;
using dim_type = std::uint8_t;
using size_type = std::size_t;

using namespace fc;

dim_type const DIM = 6;
size_type const NUM_PTS = DIM + 1;

using pc_t = point_cloud<number_type, dim_type, size_type>;
using eigen_matrix = Eigen::Matrix<number_type, Eigen::Dynamic, Eigen::Dynamic>;
using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
using eigen_cmap = Eigen::Map<const eigen_vector>;
using eigen_map = Eigen::Map<eigen_vector>;

template <typename iterator>
std::vector<number_type> project(pc_t const& pc,
                                 iterator begin, iterator end,
                                 number_type const* q) {
  std::vector<number_type> coeffs(std::distance(begin, end));
  assert(coeffs.size() > 0);
  
  eigen_matrix aff_vecs(DIM, coeffs.size() - 1);
  size_type i = 0;
  auto const orig_idx = *begin++;  // skip the origin
  while (begin != end) {
    aff_vecs.col(i) = eigen_cmap(pc[*begin], DIM) -
                      eigen_cmap(pc[orig_idx], DIM);
    i++;
    begin++;
  }
  eigen_map coeffs_map(coeffs.data(), coeffs.size());
  if (aff_vecs.cols() > 0)
    coeffs_map.tail(aff_vecs.cols()) =
    aff_vecs.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).
             solve(eigen_cmap(q, DIM));
  coeffs.front() = 1 - coeffs_map.tail(aff_vecs.cols()).sum();
  
  return std::move(coeffs);
}

int main() {
  std::srand(std::time(nullptr));
  
  eigen_matrix const data_mat = eigen_matrix::Random(NUM_PTS, DIM);
  
  std::vector<number_type const*> data_ptrs;
  for (std::size_t i = 0; i < NUM_PTS; ++i)
    data_ptrs.push_back(data_mat.row(i).data());
  pc_t pc(DIM, data_ptrs.begin(), data_ptrs.end());
                                                   
  affine_hull<number_type, size_type, dim_type> ah(pc);
  eigen_vector pt_to_project = eigen_vector::Random(DIM);
  
  // 1) add points
  std::vector<size_type> added_pts;
  eigen_vector aff_coeffs(NUM_PTS);
  for (size_type i = 0; i < pc.size(); ++i) {
    ah.add_point(i);
    added_pts.push_back(i);
    auto const coeffs = project(pc, added_pts.begin(), added_pts.end(),
                                pt_to_project.data());
    ah.project(pt_to_project.data(), aff_coeffs.data());
    assert(aff_coeffs.head(added_pts.size()).isApprox(
           eigen_cmap(coeffs.data(), added_pts.size())));
    // check the state of ah
    assert(ah.size() == added_pts.size());
    assert(std::is_permutation(ah.begin(), ah.end(), added_pts.begin()));
  }
  
  // 2) delete points, such that every sub-case is used
  while (ah.size() > 0) {
    size_type idx_to_del;
    std::vector<size_type>::iterator it_to_del;
    if (ah.size() % 2 == 1) {
      idx_to_del = added_pts.back();
      it_to_del = added_pts.begin() + (added_pts.size() - 1);
    } else {
      idx_to_del = added_pts.front();
      it_to_del = added_pts.begin();
    }
    added_pts.erase(it_to_del);
    ah.drop_point(idx_to_del);
    if (!added_pts.empty()) {
      auto const coeffs = project(pc, added_pts.begin(), added_pts.end(),
                                  pt_to_project.data());
      ah.project(pt_to_project.data(), aff_coeffs.data());
      assert(aff_coeffs.head(added_pts.size()).isApprox(
             eigen_cmap(coeffs.data(), added_pts.size())));
    }
     // check the state of ah
    assert(ah.size() == added_pts.size());
    assert(std::is_permutation(ah.begin(), ah.end(), added_pts.begin()));
  }
  
  // 3) add points after deletion
  for (size_type i = 0; i < pc.size(); ++i) {
    ah.add_point(i);
    added_pts.push_back(i);
    auto const coeffs = project(pc, added_pts.begin(), added_pts.end(),
                                pt_to_project.data());
    ah.project(pt_to_project.data(), aff_coeffs.data());
    assert(aff_coeffs.head(added_pts.size()).isApprox(
           eigen_cmap(coeffs.data(), added_pts.size())));
    // check the state of ah
    assert(ah.size() == added_pts.size());
    assert(std::is_permutation(ah.begin(), ah.end(), added_pts.begin()));
  }

  std::exit(EXIT_SUCCESS);
}
