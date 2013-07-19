#include <cstdint>

#include <iterator>
#include <type_traits>
#include <random>
#include <vector>

#include "point_cloud.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

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

TEST_CASE ("point_cloud" , "test the interface of class point_cloud") {
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

  SECTION ("ctor/asgn", "check for correct constructor/assignment behavior") {
    REQUIRE (!std::is_default_constructible<point_cloud_type>::value);
    REQUIRE (!std::is_copy_constructible<point_cloud_type>::value);
    REQUIRE (!std::is_move_constructible<point_cloud_type>::value);
    REQUIRE (!std::is_copy_assignable<point_cloud_type>::value);
    REQUIRE (!std::is_move_assignable<point_cloud_type>::value);
  }

  SECTION ("begin, end", "check for range-based for-loop compatibility") {
    REQUIRE (NUM_PTS == static_cast<size_type>(std::distance(pc.begin(), pc.end())));
    auto point_it = points.cbegin();
    for (auto && p : pc)
      REQUIRE (is_equal(p, fc_point(point_it++->data()), DIM));
  }

  SECTION ("idx-op", "correct behavior of idx-operator") {
    for (size_type i = 0; i < NUM_PTS; ++i)
      REQUIRE (is_equal(pc[i], fc_point(points[i].data()), DIM));
  }
  
  SECTION ("dim/size", "check dim and size getters") {
    REQUIRE (DIM  == pc.dim());
    REQUIRE (NUM_PTS == pc.size());
  }
}

