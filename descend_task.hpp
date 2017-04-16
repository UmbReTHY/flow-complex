#ifndef DESCEND_TASK_HPP_
#define DESCEND_TASK_HPP_

#include <ostream>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <glog/logging.h>

#include "affine_hull.hpp"
#include "common.hpp"
#include "critical_point.hpp"
#include "flow_complex.hpp"
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
  descend_task(affine_hull<point_cloud_type> ah,
               eigen_vector && location,
               critical_point<number_type, size_type> * succ)
    : _ah(std::move(ah)), _location(), _succ(succ) {
    DLOG(INFO) << "****DT-CONSTRUCTOR-BEGIN*****\n"
             << "address = " << this << std::endl
             << _ah << std::endl
             << "****DT-CONSTRUCTOR-END*****\n";
    _location.swap(location);
  }
  
  descend_task(descend_task && tmp)
    : _ah(std::move(tmp._ah)), _location(), _succ(tmp._succ),
      _ignore_indices(std::move(tmp._ignore_indices)) {
    DLOG(INFO) << "DT-MOVE-CONSTR: address = " << this << std::endl;
    _location.swap(tmp._location);
  }
  
  descend_task & operator=(descend_task && rhs) {
    DLOG(INFO) << "DT-MOVE-ASSIGN: address = " << this << std::endl;
    if (this != &rhs) {
      _ah = std::move(rhs._ah);
      _succ = rhs._succ;
      _ignore_indices = std::move(rhs._ignore_indices);
      _location.swap(rhs._location);
    }
    return *this;
  }
  
  ~descend_task() {
    DLOG(INFO) << "DELETE-DT: " << this << std::endl;
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
    eigen_vector driver(pc.dim());
    eigen_vector lambda(pc.dim() + 1);
    thread_local std::vector<size_type> idx_store(pc.size());
    std::vector<size_type> nnvec(pc.dim() + 1);
    eigen_vector ray;

    update_ray<RAY_DIR::TO_DRIVER>(_ah, _location, lambda, driver, ray);
    DLOG(INFO) << "DESCEND-TASK STARTS: address = " << this << std::endl
             << _ah << std::endl
             << "LOCATION: " << _location.transpose() << std::endl
             << "RAY: " << ray.transpose() << std::endl;
    auto nn = std::make_pair(nnvec.begin(), number_type(0));
    try {
      auto ignoreFn = [this](size_type idx) {
        return (_ah.end() != std::find(_ah.begin(), _ah.end(), idx)) ||
               (_ignore_indices.cend() !=
                std::find(_ignore_indices.cbegin(), _ignore_indices.cend(), idx));
      };
      auto vf = make_dt_filter(pc, driver, pc[*_ah.begin()], ignoreFn,
                               idx_store.begin());
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
      DLOG(INFO) << "DESCEND SUCCESSFUL\n";
      auto & x = driver;
      DCHECK(_ah.size() <= lambda.size());
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      if (std::distance(pos_offsets.begin(), pos_end) == _ah.size()) {
        DLOG(INFO) << "WITHIN CONVEX HULL\n";
        DLOG(INFO) << "LOC = " << x.transpose() << std::endl;
        number_type const sq_dist = (x - _ah.pc()[*_ah.begin()]).squaredNorm();
        cp_type new_cp(_ah.begin(), _ah.end(), sq_dist, _succ);
        const auto insert_pair =  fc.insert(std::move(new_cp));
        if (insert_pair.first) {  // it was a new critical point
          auto * new_succ = insert_pair.second;
          if (new_succ->index() > 1) {
            spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                               eigen_vector(x), _ah, new_succ);
          } else {
            // update incidences of the adjacent minima, of this gabriel edge
            DCHECK(1 == new_succ->index());
            auto idx_it = new_succ->idx_begin();
            fc.minimum(*idx_it)->add_successor(new_succ);
            fc.minimum(*++idx_it)->add_successor(new_succ);
          }
          // search for the adjacent maximum, if we found a d-1 critical point
          if (_ah.size() == pc.dim()) {
            DLOG(INFO) << "SPAWN ASCEND TO MAX ON OTHER SIDE\n";
            using at = ascend_task<point_cloud_type>;
            DLOG(INFO) << "t = " << nn.second << std::endl;
            // TODO actually we have all the information we need at this point:
            //      if (nn.first == nnvec.begin()) the maximum on the other side is
            //      inf (not true -> radius search confined search to small set of test indices ...),
            //      otherwise it's a possible max-candidate formed by all the
            //      points in _ah and *nn.first: we could just call the member
            //      simplex-case-upflow ascend task then --> saves another
            //      nn-along-ray search
            ath(at(std::move(_ah), eigen_vector(driver), std::move(ray)));
          }
        } else {
          // update incidences of insert_pair.second
          auto * already_found_cp = insert_pair.second;
          already_found_cp->add_successor(_succ);
        }
      } else {
        DLOG(INFO) << "ONLY AFFINE HULL\n";
        spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                           eigen_vector(x), _ah, _succ);
      }
    } else {
      DLOG(INFO) << "old location is " << _location.transpose();
      DCHECK(nn.second < 1.0);
      _location += nn.second * ray;
      DLOG(INFO) << "new location is " << _location.transpose();
      auto const  stopper_begin = nnvec.begin();
      auto const& stopper_end = nn.first;
      for (auto it = stopper_begin; it != stopper_end; ++it) {
        size_type stopper = *it;
        DLOG(INFO) << "DESCEND GOT STOPPED by index: " << stopper << std::endl;
        _ah.append_point(stopper);
      }
      DCHECK(lambda.size() >= _ah.size());
      // more than stopper can be < 0 -> we have to project
      _ah.project(_location, lambda.head(_ah.size()));
      DLOG(INFO) << "affine hull size = " << _ah.size();
      DLOG(INFO) << "lambda = " << lambda.head(_ah.size()).transpose();
      auto pos_end = get_pos_offsets(lambda.head(_ah.size()),
                                     pos_offsets.begin());
      spawn_sub_descends(dth, fc, pos_offsets.begin(), pos_end,
                         std::move(_location), std::move(_ah), _succ);
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

