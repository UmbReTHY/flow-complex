#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <ostream>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "affine_hull.hpp"
#include "common.hpp"
#include "critical_point.hpp"
#include <flow_complex.hpp>
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
  typedef flow_complex<number_type, size_type>   fc_type;
private:
  using eigen_vector = Eigen::Matrix<number_type, Eigen::Dynamic, 1>;
  using cp_type = critical_point<number_type, size_type>;
public:
  template <class Iterator>
  descend_task(affine_hull<point_cloud_type> ah,
               eigen_vector && location,
               critical_point<number_type, size_type> * succ,
               Iterator ignore_begin, Iterator ignore_end)
    : _ah(std::move(ah)), _location(), _succ(succ),
      _ignore_indices(ignore_begin, ignore_end) {
    Logger() << "****DT-CONSTRUCTOR-BEGIN*****\n"
             << "address = " << this << std::endl
             << _ah << std::endl
             << "****DT-CONSTRUCTOR-END*****\n";
    _location.swap(location);
  }
  
  descend_task(descend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _succ(tmp._succ),
      _ignore_indices(std::move(tmp._ignore_indices)) {
    Logger() << "DT-MOVE-CONSTR: address = " << this << std::endl;
    _location.swap(tmp._location);
  }
  
  descend_task & operator=(descend_task && rhs) {
    Logger() << "DT-MOVE-ASSIGN: address = " << this << std::endl;
    if (this != &rhs) {
      _ah = std::move(rhs._ah);
      _succ = rhs._succ;
      _ignore_indices = std::move(rhs._ignore_indices);
      _location.swap(rhs._location);
    }
    return *this;
  }
  
  ~descend_task() {
    Logger() << "DELETE-DT: " << this << std::endl;
  }
  
  descend_task(descend_task const&) = delete;
  descend_task & operator=(descend_task const&) = delete;
  
  void add_ignore_idx(size_type idx) {
    _ignore_indices.push_back(idx);
  }
  
  template <class DTHandler, class ATHandler, class CIHandler>
  void execute(ATHandler & ath, DTHandler & dth, fc_type & fc,
               CIHandler & cih) {
    auto const& pc = _ah.pc();
    thread_local eigen_vector driver(pc.dim());
    thread_local eigen_vector lambda(pc.dim() + 1);
    thread_local std::vector<size_type> idx_store(pc.size());
    thread_local std::vector<size_type> nnvec(pc.dim() + 1);
    thread_local eigen_vector ray;

    update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
    Logger() << "DESCEND-TASK STARTS: address = " << this << std::endl
             << _ah << std::endl
             << "LOCATION: " << _location.transpose() << std::endl
             << "RAY: " << ray.transpose() << std::endl;
    auto nn = std::make_pair(nnvec.begin(), number_type(0));
    try {
      auto vf = make_dt_filter(_ah, driver, _ignore_indices.begin(),
                               _ignore_indices.end(), idx_store.begin());
      size_type const max_num_nn = pc.dim() + 1 - _ah.size();
      nn = nearest_neighbor_along_ray(_location, ray, pc[*_ah.begin()], vf,
                                      nnvec.begin(),
                                      std::next(nnvec.begin(), max_num_nn));
    } catch(std::exception & e) {
      std::cerr << "DESCEND-TASK error: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    auto & pos_offsets = idx_store;
    if (nn.first == nnvec.begin() or nn.second > 1.0) {
      Logger() << "DESCEND SUCCESSFUL\n";
      auto & x = driver;
      assert(_ah.size() <= lambda.size());
      // handles both the critical and non-critical case
      using ci_type = circumsphere_ident<size_type>;
      if (_ah.size() == pc.dim() and cih(ci_type(_ah.begin(), _ah.end()))) {
        Logger() << "SPAWN ASCEND TO MAX ON OTHER SIDE\n";
        using at = ascend_task<point_cloud_type>;
        Logger() << "t = " << nn.second << std::endl;
        assert(nn.first == nnvec.begin());  // there is no stopper -> radius search
        ath(at(_ah, eigen_vector(driver), std::move(ray)));
      } else {
        Logger() << "NO ASCEND TO OTHER SIDE, BECAUSE "
                  << (_ah.size() == pc.dim() ? "ALREADY DESCENDED IN HERE\n"
                                             : "NO D-1 CRITICAL POINT\n");
      }
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      if (std::distance(pos_offsets.begin(), pos_end) == _ah.size()) {
        Logger() << "WITHIN CONVEX HULL\n";
        Logger() << "LOC = " << x.transpose() << std::endl;
        number_type const sq_dist = (x - _ah.pc()[*_ah.begin()]).squaredNorm();
        cp_type new_cp(_ah.begin(), _ah.end(), sq_dist, _succ);
        auto const insert_pair =  fc.insert(std::move(new_cp));
        if (insert_pair.first) {
          auto * new_succ = insert_pair.second;
          if (new_succ->index() > 1) {
            auto const pos_end = std::next(pos_offsets.begin(), _ah.size());
            std::iota(pos_offsets.begin(), pos_end, 0);
            spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                               eigen_vector(x), std::move(_ah), new_succ);
          } else {
            // update incidences of the adjacent minima, of this gabriel edge
            assert(1 == new_succ->index());
            auto idx_it = new_succ->idx_begin();
            fc.minimum(*idx_it)->add_successor(new_succ);
            fc.minimum(*++idx_it)->add_successor(new_succ);
          }
        } else {
          // update incidences of insert_pair.second
          auto * already_found_cp = insert_pair.second;
          already_found_cp->add_successor(_succ);
        }
      } else {
        Logger() << "ONLY AFFINE HULL\n";
        spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                           eigen_vector(x), std::move(_ah), _succ);
      }

    } else {
      assert(nn.second < 1.0);
      _location += nn.second * ray;
      auto const  stopper_begin = nnvec.begin();
      auto const& stopper_end = nn.first;
      for (auto it = std::next(stopper_begin); it != stopper_end; ++it) {
        size_type stopper = *it;
        Logger() << "DESCEND GOT STOPPED by index: " << stopper << std::endl;
        auto new_ah(_ah);
        new_ah.append_point(stopper);
        assert(lambda.size() >= new_ah.size());
        // more than stopper can be < 0 -> we have to project
        new_ah.project(_location, lambda.head(new_ah.size()));
        assert(lambda[new_ah.size() - 1] < 0);  // stopper: < 0
        auto pos_end = get_pos_offsets(lambda.head(new_ah.size()),
                                       pos_offsets.begin());
        spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                           std::move(_location), std::move(new_ah), _succ,
                           stopper_begin, stopper_end);
      }
      // move first iteration - instead of copying
      size_type const stopper = *stopper_begin;
      _ah.append_point(stopper);
      _ah.project(_location, lambda.head(_ah.size()));
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                         std::move(_location), std::move(_ah), _succ,
                         stopper_begin, stopper_end);
    }
  }
private:
  affine_hull<point_cloud_type>             _ah;
  eigen_vector                        _location;
  cp_type *                               _succ;
  std::vector<size_type>        _ignore_indices;
};

}  // namespace FC

#endif  // DESCEND_TASK_HPP_

