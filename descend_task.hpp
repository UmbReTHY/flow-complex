#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <utility>
#include <vector>

#include <Eigen/Core>

#include "critical_point.hpp"
#include "affine_hull.hpp"
#include "update_ray.hpp"
#include "vertex_filter.hpp"

namespace FC {

template <typename point_cloud_t>
class descend_task {
public:
  typedef point_cloud_t                          point_cloud_type;
  typedef typename point_cloud_type::number_type number_type;
  typedef typename point_cloud_type::size_type   size_type;
private:
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
public:
  descend_task(affine_hull<point_cloud_type> && ah,
               eigen_vector && location,
               critical_point<number_type, size_type> * succ)
    : _ah(std::move(ah)), _succ(succ) {
    _location.swap(location);
  }
  
  template <typename DTHandler, typename ATHandler, typename CPHandler>
  void execute(DTHandler & dth, ATHandler & ath, CPHandler & cph) {
    auto const& pc = _ah.pc();
    eigen_vector ray;
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
    std::vector<size_type> nnvec(pc.dim());  // for at most d additional nn
    // TODO dropped point storage does not have to be a member!
    auto vf = make_vertex_filter(_ah);
//    do {
//      auto nn = std::make_pair(nnvec.begin(), number_type(0));
//      try {
//        nn = nearest_neighbor_along_ray(_location, _ray, pc[*_ah.begin()],
//                                        vf, nnvec.begin(), nnvec.begin() +
//                                        (pc.dim() + 1 - _ah.size()));
//      } catch(std::exception & e) {
//        std::fprintf(stderr, "error: %s\n", e.what());
//        std::exit(EXIT_FAILURE);
//      }
//      using dt = descend_task<point_cloud_type>;
//      if (nn.second == std::numeric_limits<number_type>::infinity()) {
//        assert(_ah.size() == pc.dim());
//        dth(dt(std::move(_ah), std::move(_location),
//               cph(cp_type(pc.dim())).second));
//        break;  // EXIT 1
//      } else {
//        _location += nn.second * _ray;
//        for (auto it = nnvec.begin(); it != nn.first; ++it)
//          _ah.add_point(*it);
//        // check for finite max
//        if (_ah.size() == pc.dim() + 1) {
//          eigen_vector lambda(_ah.size());
//          _ah.project(_location, lambda);
//          // drop negative indices
//          auto m_it = _ah.begin();
//          for (size_type i = 0; i < lambda.size(); ++i) {
//            // TODO remember dropped indices
//            if (lambda[i] < 0)
//              _ah.drop_point(*m_it);
//            m_it++;
//          }
//          if (_ah.size() == pc.dim() + 1) {
//            assert(_ah.size() == pc.dim() + 1);
//            number_type dist = (_location - pc[*_ah.begin()]).norm();
//            auto r_pair = cph(cp_type(_ah.begin(), _ah.end(), dist));
//            if (r_pair.first) {
//              auto max_ptr = r_pair.second;
//              if (max_ptr->index() > 1) { // don't descend to idx-0 cp-s
//                // TODO maybe pass the recently dropped indices
//                for (auto it = ++_ah.begin(); it != _ah.end(); ++it) {
//                  auto new_ah = _ah;
//                  new_ah.drop_point(*it);
//                  dth(dt(std::move(new_ah), eigen_vector(_location), max_ptr));
//                }
//                // reuse this task's affine hull for one of the descent tasks
//                _ah.drop_point(*_ah.begin());
//                dth(dt(std::move(_ah), std::move(_location), max_ptr));
//              }
//            }
//            break;    // EXIT 2 - none have been dropped -> finite max
//          } else {
//            update_ray(_ah, _location, lambda, driver, _ray);
//          }
//        } else {
//          update_ray(_ah, _location, lambda, driver, _ray);
//        }
//      }
//      vf.reset();  // make the vertex filter consider all points
//                   // of the point cloud again on subsequent calls
//    } while(true);
  }
  
private:
  affine_hull<point_cloud_type>            _ah;
  eigen_vector                             _location;
  critical_point<number_type, size_type> * _succ;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

