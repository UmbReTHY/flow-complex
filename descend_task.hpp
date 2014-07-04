#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <utility>
#include <vector>

#include <Eigen/Core>

#include "common.hpp"
#include "critical_point.hpp"
#include "affine_hull.hpp"
#include "update_ray.hpp"
#include "vertex_filter.hpp"
#include "utility.hpp"

namespace FC {

template <typename point_cloud_t>
class ascend_task;  // break cyclic dependencies

template <class Derived, class Iterator>  // forward declaration
Iterator get_pos_offsets(Eigen::MatrixBase<Derived> const& , Iterator );

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
  descend_task(affine_hull<point_cloud_type> ah,
               eigen_vector && location,
               critical_point<number_type, size_type> * succ,
               size_type ignore_idx)
    : _ah(std::move(ah)), _location(), _succ(succ), _ignore_idx(ignore_idx) {

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
      _ignore_idx(tmp._ignore_idx) {

    std::cout << "DT-MOVE-CONSTR: address = " << this << std::endl;

    _location.swap(tmp._location);
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
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    std::vector<size_type> idx_store(pc.dim() + 1);
    eigen_vector ray;
    update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
    
    static int dt_id = 0;
    dt_id++;
    std::cout << "DESCEND_TASK_ID = " << dt_id << std::endl;
    std::cout << "address = " << this << std::endl;
    std::cout << _ah << std::endl;
    std::cout << "LOCATION: " << _location.transpose() << std::endl;
    std::cout << "RAY: " << ray.transpose() << std::endl;

    size_type stopper;  // dts asumme at most 1 stopper: non-degenerate case!
    auto nn = std::make_pair(&stopper, number_type(0));
    try {
      auto vf = make_vertex_filter(_ah, &_ignore_idx, std::next(&_ignore_idx));
      nn = nearest_neighbor_along_ray(_location, ray, pc[*_ah.begin()], vf,
                                      &stopper, std::next(&stopper));
    } catch(std::exception & e) {
      std::fprintf(stderr, "error: %s\n", e.what());
      std::exit(EXIT_FAILURE);
    }
    auto & pos_offsets = idx_store;
    if (nn.first == &stopper or nn.second > 1.0) {
      std::cout << "DESCEND SUCCESSFUL\n";
      auto & x = driver;
      assert(_ah.size() <= lambda.size());
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      if (std::distance(pos_offsets.begin(), pos_end) == _ah.size()) {
        std::cout << "WITHIN CONVEX HULL\n";
        std::cout << "LOC = " << x.transpose() << std::endl;
        number_type const sq_dist = (x - _ah.pc()[*_ah.begin()]).squaredNorm();
        auto const insert_pair =  cph(cp_type(_ah.begin(), _ah.end(), sq_dist,
                                      _succ));
        if (insert_pair.first) {
          auto * new_succ = insert_pair.second;
          if (new_succ->index() > 1) {
            // TODO for _ah.size() != pc.dim(), we can move the arguments
            std::iota(pos_offsets.begin(),
                      std::next(pos_offsets.begin(), _ah.size()), 0);
            spawn_sub_descends(dth, pos_offsets.begin(),
                               std::next(pos_offsets.begin(), _ah.size()),
                               eigen_vector(x), _ah, new_succ);
          } else {
            // TODO update incidences of idx 0 cps
          }
        } else {
          // TODO update incidences of insert_pair.second
        }
      } else {
        std::cout << "ONLY AFFINE HULL\n";
        // TODO move _ah and x if the ASCEND branch below won't be executed
        spawn_sub_descends(dth, pos_offsets.begin(), pos_end,
                           eigen_vector(x), _ah, _succ);
      }
      // handles both the critical and non-critical case
      using ci_type = circumsphere_ident<size_type>;
      if (_ah.size() == pc.dim() and cih(ci_type(_ah.begin(), _ah.end()))) {
        std::cout << "SPAWN ASCEND TO MAX ON OTHER SIDE\n";
        using at = ascend_task<point_cloud_type>;
        auto const& t = nn.second;
        if (nn.first != &stopper) { // there is a stopper
          assert(t > 1.0);
          _location += t * ray;
          _ah.add_point(stopper);
          // TODO if this routine could return weather the simplex was critical or not,
          //      we could add just another succ to spawned sub descends spawned from here
          //      an avoid re-traversal of descends back to where this new maxima came from
          simplex_case_upflow(std::move(_ah), eigen_vector(_location), lambda,
                              std::move(ray), driver, idx_store.begin(),
                              idx_store.end(), cph, ath, dth);
        } else {
          // TODO introduce vector of succs -> for alle spawn_sub_descends
          //      we just have to update this tasks succ vector
          //      (with the cp at inf) and pass it; we have to move this check
          //      way up to the top, however
        }
      } else {
        std::cout << "NO ASCEND TO OTHER SIDE, BECAUSE "
                  << (_ah.size() == pc.dim() ? "ALREADY DESCENDED IN HERE\n"
                                             : "NO D-1 CRITICAL POINT\n");
      }
    } else {
      std::cout << "DESCEND GOT STOPPED by index " << stopper << std::endl;
      assert(nn.second < 1.0);
      _location += nn.second * ray;
      _ah.add_point(stopper);  // TODO rename to append_point
      assert(lambda.size() >= _ah.size());
      _ah.project(_location, lambda.head(_ah.size()));  // more than stopper can be < 0
      assert(lambda[_ah.size() - 1] < 0);  // stopper: < 0
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      spawn_sub_descends(dth, pos_offsets.begin(), pos_end,
                         std::move(_location), std::move(_ah), _succ);
    }
  }
  
private:
  ah_type       _ah;
  eigen_vector  _location;
  cp_type *     _succ;
  size_type     _ignore_idx;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

