#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <utility>
#include <vector>

#include <Eigen/Core>

#include "critical_point.hpp"
#include "affine_hull.hpp"
#include "update_ray.hpp"
#include "vertex_filter.hpp"
#include "utility.hpp"

namespace FC {

template <typename point_cloud_t>
class ascend_task;  // break cyclic dependencies

template <typename point_cloud_t>
class descend_task {
public:
  typedef point_cloud_t                          point_cloud_type;
  typedef typename point_cloud_type::number_type number_type;
  typedef typename point_cloud_type::size_type   size_type;
private:
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using cp_type = critical_point<number_type, size_type>;
  using ah_type = affine_hull<point_cloud_type>;
public:
  descend_task(affine_hull<point_cloud_type> && ah,
               eigen_vector && location,
               critical_point<number_type, size_type> * succ)
    : _ah(std::move(ah)), _succ(succ) {
    std::cout << "****DT-CONSTRUCTOR-BEGIN*****\n";
    std::cout << "address = " << this << std::endl;
    std::cout << "HULL MEMBERS: ";
    for (auto it = _ah.begin(); it != _ah.end(); ++it)
      std::cout << *it << ", ";
    std::cout << std::endl;
    std::cout << "****DT-CONSTRUCTOR-END*****\n";
    _location.swap(location);
  }
  
  descend_task(descend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _succ(tmp._succ) {
    _location.swap(tmp._location);
    std::cout << "DT-MOVE-CONSTR: address = " << this << std::endl;
  }
  
  ~descend_task() {
    std::cout << "DELETE-DT: " << this << std::endl;
  }
  
  descend_task(descend_task const&) = delete;
  descend_task & operator=(descend_task const&) = delete;
  
  template <typename DTHandler, typename ATHandler, typename CPHandler>
  void execute(DTHandler & dth, ATHandler & ath, CPHandler & cph) {
    auto const& pc = _ah.pc();
    // TODO the vectors below would all qualify for thread local storage
    eigen_vector ray;
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    std::vector<size_type> nnvec(pc.dim());  // for at most d additional nn
    std::vector<size_type> dropvec(pc.dim());
    auto dropped_end = dropvec.begin();
    using idx_iterator = typename affine_hull<point_cloud_type>::iterator;
    auto vf = make_vertex_filter(_ah, dropvec.begin(), dropped_end);
    
    static int dt_id = 0;
    dt_id++;
    std::cout << "DESCEND_TASK_ID = " << dt_id << std::endl;
    std::cout << "address = " << this << std::endl;
    
    do {
      std::cout << "_ah.size() = " << _ah.size() << std::endl;
      std::cout << "HULL MEMBERS: ";
      for (auto it = _ah.begin(); it != _ah.end(); ++it)
        std::cout << *it << ", ";
      std::cout << std::endl;
      std::cout << "LOCATION: " << _location.transpose() << std::endl;
    
      update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
      std::cout << "RAY: " << ray.transpose() << std::endl;
      
      auto nn = std::make_pair(nnvec.begin(), number_type(0));
      try {
        std::cout << "NORM OF RAY " << ray.squaredNorm() << std::endl;
        nn = nearest_neighbor_along_ray(_location, ray, pc[*_ah.begin()],
                                        vf, nnvec.begin(), nnvec.begin() +
                                        (pc.dim() + 1 - _ah.size()));
                                        
//        std::cout << "t = " << nn.second << std::endl;
                                        
      } catch(std::exception & e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        std::exit(EXIT_FAILURE);
      }
      if (nn.first == nnvec.begin() or nn.second > 1.0) {
        std::cout << "DESCEND SUCCESSFUL\n";
        _location = driver;  // TODO this assignment can be deferred until after the next if
        // drop neg coeffs
        dropped_end = drop_neg_coeffs(_location, lambda.head(_ah.size()), _ah,
                                      dropvec.begin());  // TODO we know the negative coeffs already
                                                         //      from the computation of the driver
        if (dropped_end == dropvec.begin()) {
          std::cout << "WITHIN CONVEX HULL\n";
          std::cout << "LOC = " << _location.transpose() << std::endl;
          auto const insert_pair = 
          cph(cp_type(_ah.begin(), _ah.end(),
                      (_location - _ah.pc()[*_ah.begin()]).squaredNorm(),
                      _succ));
          if (insert_pair.first) {
            if (insert_pair.second->index() == (pc.dim() - 1)) {
              using at = ascend_task<point_cloud_type>;
              ath(at(_ah, _location, ray));  // TODO for dim() == 2, we can move the arguments
            }
            if (insert_pair.second->index() > 2) {
              spawn_sub_descends(dth, _location, _ah.size(),
                                 insert_pair.second, _ah);
              _succ = insert_pair.second;  // for all subsequent discoveries
                                           // we have a new direct succ

              std::cout << "SPAWNING(dt) DT with MEMBERS: ";   
              for (auto it = _ah.begin(); it != _ah.end(); ++it)
                std::cout << *it << ", ";
              std::cout << std::endl;  // TODO we could remember the dropped
                                       //      vertex if spwan_sub_descends would not alter _ah
            } else {
              break;  // EXIT 1
            }
          } else {
            // TODO test if it is necessary for correctness to ascend in the
            //      d - 1 cp even though the cp was already found
            break;  // EXIT 2
          }
        } else {
          std::cout << "WITHIN AFFINE HULL\n";
        }
      } else {
        std::cout << "DESCEND GOT STOPPED\n";  // TODO: CONTINUE TRACKING HERE
        assert(nn.second < 1.0);
        _location += nn.second * ray;  // TODO: since norm(ray) == dist-to-driver, we can remove the multiplication here
        size_type const size_before = _ah.size();
        std::cout << "_ah.size() = " << _ah.size() << std::endl;
        for (auto it = nnvec.begin(); it != nn.first; ++it) {
          _ah.add_point(*it);  // TODO convert add_point to append_point,
                               //      to fix that semantic. This influences
                               //      correctness of spawn_sub_descends
          std::cout << "STOPPER-ID = " << *it << std::endl;
        }
        std::cout << "_ah.size() = " << _ah.size() << std::endl;
        spawn_sub_descends(dth, _location, size_before, _succ, _ah);
        std::cout << "after sub-spawn: _ah.size() = " << _ah.size() << std::endl;
      }
      vf.reset(dropped_end);  // make the vertex filter consider all points
                              // of the point cloud again on subsequent calls
      std::cout << "ONE MORE TURN\n";
      std::cout << "SPAWNING(dt) DT with MEMBERS: ";   
      for (auto it = _ah.begin(); it != _ah.end(); ++it)
        std::cout << *it << ", ";
      std::cout << std::endl;
    } while(true);
  }
  
private:
  ah_type      _ah;
  eigen_vector _location;
  cp_type *    _succ;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

