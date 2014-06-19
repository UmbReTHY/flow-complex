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
  using idx_cont = std::vector<size_type>;
public:
  template <class Iterator>
  descend_task(affine_hull<point_cloud_type> ah,
               eigen_vector location,
               critical_point<number_type, size_type> * succ,
               Iterator drop_begin, Iterator drop_end)
    : _ah(std::move(ah)), _succ(succ), _dropvec(drop_begin, drop_end) {
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
    : _ah(std::move(tmp._ah)), _location(), _succ(tmp._succ),
      _dropvec(std::move(tmp._dropvec)) {
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
    typename idx_cont::iterator dropped_end;
    { auto offset = _dropvec.size();
      _dropvec.resize(pc.dim());
      dropped_end = _dropvec.begin() + offset;
    }    
    using idx_iterator = typename affine_hull<point_cloud_type>::iterator;
    auto vf = make_vertex_filter(_ah, _dropvec.begin(), dropped_end);
    
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
      } catch(std::exception & e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        std::exit(EXIT_FAILURE);
      }
      if (nn.first == nnvec.begin() or nn.second > 1.0) {
        std::cout << "DESCEND SUCCESSFUL\n";
        _location = driver;
        dropped_end = drop_neg_coeffs(_location, lambda.head(_ah.size()), _ah,
                                      _dropvec.begin());  // TODO we know the negative coeffs already
                                                          //      from the computation of the driver
        if (dropped_end == _dropvec.begin()) {
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
            if (insert_pair.second->index() > 1)
              spawn_sub_descends(dth, _ah.size(), std::move(_location),
                                 std::move(_ah), insert_pair.second);
          } else {
            // TODO test if it is necessary for correctness to ascend in the
            //      d - 1 cp even though the cp was already found
          }
          break;  // EXIT 1
        } else {
          std::cout << "WITHIN AFFINE HULL - ONE MORE TURN\n";
          vf.reset(dropped_end);  // make the vertex filter consider all points
                                  // of the point cloud again on subsequent calls
          std::cout << "SPAWNING(dt) DT with MEMBERS: ";   
          for (auto it = _ah.begin(); it != _ah.end(); ++it)
            std::cout << *it << ", ";
          std::cout << std::endl;
        }
      } else {
        std::cout << "DESCEND GOT STOPPED\n";  // TODO: CONTINUE TRACKING HERE
        assert(nn.second < 1.0);
        _location += nn.second * ray;
        size_type const size_before = _ah.size();
        std::cout << "_ah.size() = " << _ah.size() << std::endl;
        for (auto it = nnvec.begin(); it != nn.first; ++it) {
          _ah.add_point(*it);  // TODO convert add_point to append_point,
                               //      to fix that semantic. This influences
                               //      correctness of spawn_sub_descends
          std::cout << "STOPPER-ID = " << *it << std::endl;
        }
        spawn_sub_descends(dth, size_before, std::move(_location),
                           std::move(_ah), _succ);
        break;  // EXIT 2
      }
    } while(true);
  }
  
private:
  ah_type       _ah;
  eigen_vector  _location;
  cp_type *     _succ;
  idx_cont      _dropvec;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

