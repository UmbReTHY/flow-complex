#ifndef ASCEND_TASK_HPP_
#define ASCEND_TASK_HPP_

#include <utility>
#include <random>
#include <tuple>

#include <Eigen/Dense>

#include "descend_task.hpp"
#include "point_cloud.hpp"
#include "affine_hull.hpp"

namespace FC {

template <typename _number_type, typename _size_type>
class ascend_task {
  using eigen_vector = Eigen::Matrix<_number_type, Eigen::Dynamic, 1>;

  public:
    typedef _number_type                             number_type;
    typedef _size_type                                 size_type;
    typedef point_cloud<number_type, size_type> point_cloud_type;
    typedef descend_task<number_type, size_type> descend_task_type;
    
    ascend_task(point_cloud_type const& pc)
      : _nn_aff_hull(pc), _location(pc.dim()), _ray(pc.dim()) {
      using cmap = Eigen::Map<eigen_vector const>;
      // generate seed
      std::tuple<_size_type, _number_type, bool> nn;
      // reseed as long as the nearest neighbor is not unique
      do {
        gen_convex_comb(pc, _location);
        nn = pc.nearest_neighbor(&_location);
      } while (std::get<2>(nn));
      // add the nearest neighbor
      _nn_aff_hull.add_point(std::get<0>(nn));
      // set the ray
      _ray = cmap(pc[*_nn_aff_hull.begin()], pc.dim()) - _location;
    }
    
    ascend_task(descend_task_type && dt)
    : _nn_aff_hull(std::move(dt._nn_aff_hull)) {
      _location.swap(dt._location);
      _ray.swap(dt._ray);
    }
    
    /**
      @return true, if finite maximum, false if maximum at infinity
    */
    bool execute() {
      // TODO continue
      return true;
    }
    
  private:
    /**
      @brief generates a convex combination of a set of points, where the
             coefficients of the convex combination are strictly greater than 0
      @param begin, end iterators to pointers that hold the coeeficients of the
             points
      @param dim the the size of the point arrays
      @param target a pointer to at least dim elements where the result is placed
             
    */
    void gen_convex_comb(point_cloud_type const& pc, eigen_vector & target) {
      std::random_device rd;
      std::mt19937 gen(rd());
      // it's important for the lower bound not to be 0, otherwise some points
      // might not take part in the convex combination and the seed point would be
      // located on the convex hull, thus getting stuck during subsequent ascends
      // and not reaching a maximum
      std::uniform_real_distribution<number_type> dis(0.5, 1.0);
      number_type sum = 0.0;
      target.setZero();
      number_type tmp;
      for (auto it = pc.cbegin(); it != pc.cend(); ++it) {
        tmp = dis(gen);
        sum += tmp;
        target = tmp * Eigen::Map<eigen_vector const>(*it, target.size());
      }
      target /= sum;
    }
  
    affine_hull<_number_type, _size_type> _nn_aff_hull;
    eigen_vector                          _location;
    eigen_vector                          _ray;
    // TODO we probably need a store for dropped points
};

}  // namespace

#endif  // ASCEND_TASK_HPP_

