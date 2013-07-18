#include <cassert>
#include <cstdlib>
#include <cstdint>

#include <iterator>
#include <type_traits>
#include <random>
#include <vector>

#include <catch.hpp>

#include "point_cloud.hpp"

using number_type = double;
using dim_type = std::uint8_t;
using size_type = std::uint8_t;
using my_point = std::vector<number_type>;
using my_point_set = std::vector<my_point>;
using fc_point = FC::point<number_type>;

bool is_equal(fc_point const& a, fc_point const& b,
              dim_type dim) {
  for (dim_type i = 0; i < dim; ++i)
    if (a[i] != b[i])
      return false;
  return true;
}

int main(int, char**) {
  // initialize random distributions
  std::random_device rd;  // provides the seed
  std::mt19937 gen(rd());  // random number generator
  std::uniform_real_distribution<number_type> real_dis(0, 1);
  std::uniform_int_distribution<dim_type> dim_dis(1, 20);
  std::uniform_int_distribution<size_type> num_pts_dis(50, 100);

  // initialize external points
  dim_type const DIM = dim_dis(gen);
  size_type const NUM_PTS = num_pts_dis(gen);


  // create and initialize and the points
  my_point_set points(NUM_PTS, my_point(DIM));
  for (auto & p : points)
    for (auto & el : p)
      el = real_dis(gen);

  // construct point cloud adaptor
  using point_cloud_type = FC::point_cloud<number_type, dim_type, size_type>;
  point_cloud_type pc(DIM, points.cbegin(), points.cend());

  // 1) check for correct constructor/assignment behavior
  assert(!std::is_default_constructible<point_cloud_type>::value);
  assert(!std::is_copy_constructible<point_cloud_type>::value);
  assert(!std::is_move_constructible<point_cloud_type>::value);
  assert(!std::is_copy_assignable<point_cloud_type>::value);
  assert(!std::is_move_assignable<point_cloud_type>::value);

  // 2) check begin, end
  assert(NUM_PTS == static_cast<size_type>(std::distance(pc.begin(), pc.end())));
  auto point_it = points.cbegin();
  for (auto && p : pc)
    assert(is_equal(p, fc_point(point_it++->data()), DIM));

  // 3) check indexing operator
  for (size_type i = 0; i < NUM_PTS; ++i)
    assert(is_equal(pc[i], fc_point(points[i].data()), DIM));

	// 4) check dim() and size()
  assert(DIM  == pc.dim());
  assert(NUM_PTS == pc.size());

  std::exit(EXIT_SUCCESS);
}
