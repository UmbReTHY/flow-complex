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

template <typename _size_type>
class circumsphere_ident;  // break cyclic dependencies

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
               eigen_vector && location,
               critical_point<number_type, size_type> * succ,
               Iterator drop_begin, Iterator drop_end)
    : _ah(std::move(ah)), _location(), _succ(succ),
      _dropvec(drop_begin, drop_end) {
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
  
  template <class DTHandler, class ATHandler, class CPHandler, class CIHandler>
  void execute(DTHandler & dth, ATHandler & ath, CPHandler & cph,
               CIHandler & cih) {
    auto const& pc = _ah.pc();
    // TODO the vectors below would all qualify for thread local storage
    eigen_vector ray;
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
    
    static int dt_id = 0;
    dt_id++;
    std::cout << "DESCEND_TASK_ID = " << dt_id << std::endl;
    std::cout << "address = " << this << std::endl;
    std::cout << "HULL MEMBERS: ";
    for (auto it = _ah.begin(); it != _ah.end(); ++it)
      std::cout << *it << ", ";
    std::cout << std::endl;
    std::cout << "LOCATION: " << _location.transpose() << std::endl;
    std::cout << "RAY: " << ray.transpose() << std::endl;

    size_type stopper;  // dts asumme at most 1 stopper: non-degenerate case!    
    auto nn = std::make_pair(&stopper, number_type(0));
    try {
      auto vf = make_vertex_filter(_ah, _dropvec.begin(), _dropvec.end());
      nn = nearest_neighbor_along_ray(_location, ray, pc[*_ah.begin()],
                                      vf, &stopper, &stopper + 1);
    } catch(std::exception & e) {
      std::fprintf(stderr, "error: %s\n", e.what());
      std::exit(EXIT_FAILURE);
    }
    auto & pos_offsets = _dropvec;  // reuse _dropvec
    if (nn.first == &stopper or nn.second > 1.0) {
      std::cout << "DESCEND SUCCESSFUL\n";
      using ci_type = circumsphere_ident<size_type>;
      if (cih(ci_type(_ah.begin(), _ah.end()))) {
        auto & x = driver;
        if (_ah.size() == pc.dim()) {  // handles the dangling case as well
          std::cout << "SPAWN ASCEND\n";
          using at = ascend_task<point_cloud_type>;
          ath(at(_ah, x, ray));  // TODO for dim() == 2, we can move the arguments
                                 // also: we know the nn-along =-ray in the next
                                 // ascend already, in nn.first  != &stopper;
                                 // think of a way to pass that information,
                                 // or better: start execute at the right subcase directly
                                 // TODO: this can ping pong: we need a store for
                                 //       the non-critical (d-1, for now) affine hulls
        }
        assert(_ah.size() <= lambda.size());
        // erase positions of negative coefficients
        pos_offsets.clear();
        pos_offsets.reserve(_ah.size());
        for (size_type i = 0; i < _ah.size(); ++i) {
          if (lambda[i] > 0) { // TODO externalize
            pos_offsets.push_back(i);
            std::cout << "POS-COEFF " << lambda[i] << " of index " << *(_ah.begin() + i) << std::endl;
          } else {
            std::cout << "NEG-COEFF " << lambda[i] << " of index " << *(_ah.begin() + i) << std::endl;
          }
        }
        if (pos_offsets.size() == _ah.size()) {
          std::cout << "WITHIN CONVEX HULL\n";
          std::cout << "LOC = " << x.transpose() << std::endl;
          auto const insert_pair = 
          cph(cp_type(_ah.begin(), _ah.end(),
                      (x - _ah.pc()[*_ah.begin()]).squaredNorm(), _succ));
          if (insert_pair.first) {
            auto * new_succ = insert_pair.second;
            if (new_succ->index() > 1) {
              spawn_sub_descends(dth, pos_offsets.begin(), pos_offsets.end(),
                                 std::move(x), std::move(_ah), new_succ);
            } else {
              // TODO update incidences of idx 0 cps
            }
          }
        } else {
          std::cout << "ONLY AFFINE HULL\n";
          spawn_sub_descends(dth, pos_offsets.begin(), pos_offsets.end(),
                             std::move(x), std::move(_ah), _succ);
        }
      }
    } else {
      std::cout << "DESCEND GOT STOPPED by index " << stopper << std::endl;
      assert(nn.second < 1.0);
      _location += nn.second * ray;
//      std::iota(pos_offsets.begin(), pos_offsets.end(), 0);  // TODO WRONG: not only the stopper can be negative
      pos_offsets.clear();  // TODO duplicate code (partly)
      pos_offsets.reserve(_ah.size());
      _ah.add_point(stopper);  // TODO rename to append_point
      assert(lambda.size() >= _ah.size());
      _ah.project(_location, lambda.head(_ah.size()));  // more than stopper can be < 0
      assert(lambda[_ah.size() - 1] < 0);  // stopper: < 0
      for (size_type i = 0; i < _ah.size(); ++i) {
        if (lambda[i] > 0) { // TODO externalize
          pos_offsets.push_back(i);
          std::cout << "POS-COEFF " << lambda[i] << " of index " << *(_ah.begin() + i) << std::endl;
        } else {
          std::cout << "NEG-COEFF " << lambda[i] << " of index " << *(_ah.begin() + i) << std::endl;
        }
      }
      spawn_sub_descends(dth, pos_offsets.begin(), pos_offsets.end(),
                         std::move(_location), std::move(_ah), _succ);
    }
  }
  
private:
  ah_type       _ah;
  eigen_vector  _location;
  cp_type *     _succ;
  idx_cont     _dropvec;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

