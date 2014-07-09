#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <ostream>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "affine_hull.hpp"
#include "common.hpp"
#include "critical_point.hpp"
#include "update_ray.hpp"
#include "utility.hpp"
#include "vertex_filter.hpp"

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
public:

  descend_task(affine_hull<point_cloud_type> ah,
               eigen_vector && location,
               critical_point<number_type, size_type> const* succ,
               size_type ignore_idx)
    : _ah(std::move(ah)), _location(), _succ(succ), _ignore_idx(ignore_idx) {
    Logger() << "****DT-CONSTRUCTOR-BEGIN*****\n"
             << "address = " << this << std::endl
             << _ah << std::endl
             << "****DT-CONSTRUCTOR-END*****\n";
    _location.swap(location);
  }
  
  descend_task(descend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _succ(tmp._succ),
      _ignore_idx(tmp._ignore_idx) {
    Logger() << "DT-MOVE-CONSTR: address = " << this << std::endl;
    _location.swap(tmp._location);
  }
  
  ~descend_task() {
    Logger() << "DELETE-DT: " << this << std::endl;
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
    // TODO maybe collaps with idx_store
    std::vector<size_type> bigger_store(pc.size());
    eigen_vector ray;
    update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
    
    Logger() << "DESCEND-TASK STARTS: address = " << this << std::endl
             << _ah << std::endl
             << "LOCATION: " << _location.transpose() << std::endl
             << "RAY: " << ray.transpose() << std::endl;
    size_type stopper;  // dts asumme at most 1 stopper: non-degenerate case!
    auto nn = std::make_pair(&stopper, number_type(0));
    try {
      auto vf = make_vertex_filter(_ah, driver, _ignore_idx, bigger_store.begin());
      nn = nearest_neighbor_along_ray(_location, ray, pc[*_ah.begin()], vf,
                                      &stopper, std::next(&stopper));
    } catch(std::exception & e) {
      std::cerr << "error: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto & pos_offsets = idx_store;
    if (nn.first == &stopper or nn.second > 1.0) {
      Logger() << "DESCEND SUCCESSFUL\n";
      auto & x = driver;
      assert(_ah.size() <= lambda.size());
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      if (std::distance(pos_offsets.begin(), pos_end) == _ah.size()) {
        Logger() << "WITHIN CONVEX HULL\n";
        Logger() << "LOC = " << x.transpose() << std::endl;
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
        Logger() << "ONLY AFFINE HULL\n";
        // TODO move _ah and x if the ASCEND branch below won't be executed
        spawn_sub_descends(dth, pos_offsets.begin(), pos_end,
                           eigen_vector(x), _ah, _succ);
      }
      // handles both the critical and non-critical case
      using ci_type = circumsphere_ident<size_type>;
      if (_ah.size() == pc.dim() and cih(ci_type(_ah.begin(), _ah.end()))) {
        Logger() << "SPAWN ASCEND TO MAX ON OTHER SIDE\n";
        using at = ascend_task<point_cloud_type>;
        auto const& t = nn.second;
        assert(nn.first == &stopper);  // there is no stopper -> radius search
        ath(at(std::move(_ah), std::move(driver), std::move(ray)));
        // TODO this case has changed, especially when the subsequend at flows to inf
        // TODO introduce vector of succs -> for alle spawn_sub_descends
        //      we just have to update this tasks succ vector
        //      (with the cp at inf) and pass it; we have to move this check
        //      way up to the top, however
      } else {
        Logger() << "NO ASCEND TO OTHER SIDE, BECAUSE "
                  << (_ah.size() == pc.dim() ? "ALREADY DESCENDED IN HERE\n"
                                             : "NO D-1 CRITICAL POINT\n");
      }
    } else {
      Logger() << "DESCEND GOT STOPPED by index " << stopper << std::endl;
      assert(nn.second < 1.0);
      _location += nn.second * ray;
      _ah.add_point(stopper);  // TODO rename to append_point
      assert(lambda.size() >= _ah.size());
      // more than stopper can be < 0 -> we have to project
      _ah.project(_location, lambda.head(_ah.size()));
      assert(lambda[_ah.size() - 1] < 0);  // stopper: < 0
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      spawn_sub_descends(dth, pos_offsets.begin(), pos_end,
                         std::move(_location), std::move(_ah), _succ);
    }
  }
  
private:
  affine_hull<point_cloud_type>         _ah;
  eigen_vector                    _location;
  cp_type const*                      _succ;
  size_type                     _ignore_idx;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

