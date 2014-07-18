#ifndef TBB_LL_HPP_
#define TBB_LL_HPP_

#include <functional>
#include <utility>
#include <vector>

#include <tbb/task.h>

#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "utility.hpp"

namespace FC {

template <class pc_type_> class tbb_dt_task;

template <class pc_type_>
class tbb_at_task : public tbb::task {
public:
  typedef pc_type_                                 pc_type;
  typedef ascend_task<pc_type>                     at_type;
  typedef descend_task<pc_type>                    dt_type;
  typedef typename pc_type::number_type        number_type;
  typedef typename pc_type::size_type            size_type;
  typedef flow_complex<number_type, size_type>     fc_type;
  typedef circumsphere_ident<size_type>            ci_type;
  typedef std::function<bool(ci_type)>            cih_type;

  tbb_at_task(at_type at, fc_type & fc, cih_type cih)
  : _at(std::move(at)), _fc(fc), _cih(std::move(cih)) {
  }

  virtual tbb::task * execute() {
    size_type const dim = _fc.max_at_inf()->index();
    std::vector<at_type> at_tasks;
    at_tasks.reserve(dim);
    std::vector<dt_type> dt_tasks;
    dt_tasks.reserve(dim + 1);
    // lambdas for piling up tasks
    auto ath = [&](at_type && at) {at_tasks.push_back(std::move(at));};
    auto dth = [&](dt_type && dt) {dt_tasks.push_back(std::move(dt));};
    // perform actual task
    _at.execute(ath, dth, _fc, _cih);
    // recycling and spawning
    tbb::task_list at_list;
    tbb::task * next_task = nullptr;
    tbb::task * p = this;
    if (not at_tasks.empty()) {
      p = new(allocate_continuation())tbb::empty_task();
      p->set_ref_count(at_tasks.size() + dt_tasks.size());
      for (auto it = std::next(at_tasks.begin()); it != at_tasks.end(); ++it) {
        using tbb_at_t = tbb_at_task<pc_type>;
        auto * t = new(p->allocate_child()) tbb_at_t(std::move(*it), _fc, _cih);
        at_list.push_back(*t);
      }
      spawn(at_list);
      recycle_as_child_of(*p);
      _at = std::move(at_tasks.front());
      next_task = this;
    }
    tbb::task_list dt_list;
    for (auto && dt : dt_tasks) {
      using tbb_dt_t = tbb_dt_task<pc_type>; 
      auto * t = new(p->allocate_child()) tbb_dt_t(std::move(dt), _fc, _cih);
      dt_list.push_back(*t);
    }
    if (at_tasks.empty()) {  // we couldn't recycle anything
      set_ref_count(dt_tasks.size() + 1);
      spawn_and_wait_for_all(dt_list);
    } else {
      spawn(dt_list);
    }
    return next_task;
  }
private:
  at_type    _at;
  fc_type &  _fc;
  cih_type  _cih;
};

template <class pc_type_>
class tbb_dt_task : public tbb::task {
public:
  typedef pc_type_                                 pc_type;
  typedef ascend_task<pc_type>                     at_type;
  typedef descend_task<pc_type>                    dt_type;
  typedef typename pc_type::number_type        number_type;
  typedef typename pc_type::size_type            size_type;
  typedef flow_complex<number_type, size_type>     fc_type;
  typedef circumsphere_ident<size_type>            ci_type;
  typedef std::function<bool(ci_type)>            cih_type;

  tbb_dt_task(dt_type dt, fc_type & fc, cih_type cih)
  : _dt(std::move(dt)), _fc(fc), _cih(std::move(cih)) {
  }

  virtual tbb::task * execute() {
    size_type const dim = _fc.max_at_inf()->index();
    std::vector<at_type> at_tasks;
    at_tasks.reserve(dim);
    std::vector<dt_type> dt_tasks;
    dt_tasks.reserve(dim + 1);
    // lambdas for piling up tasks
    auto ath = [&](at_type && at) {at_tasks.push_back(std::move(at));};
    auto dth = [&](dt_type && dt) {dt_tasks.push_back(std::move(dt));};
    // perform actual task
    _dt.execute(ath, dth, _fc, _cih);
    // recycling and spawning
    tbb::task * next_task = nullptr;
    tbb::task * p = this;
    tbb::task_list dt_list;
    if (not dt_tasks.empty()) {
      p = new(allocate_continuation())tbb::empty_task();
      p->set_ref_count(at_tasks.size() + dt_tasks.size());
      for (auto it = std::next(dt_tasks.begin()); it != dt_tasks.end(); ++it) {
        using tbb_dt_t = tbb_dt_task<pc_type>;
        auto * t = new(p->allocate_child()) tbb_dt_t(std::move(*it), _fc, _cih);
        dt_list.push_back(*t);
      }
      spawn(dt_list);
      recycle_as_child_of(*p);
      _dt = std::move(dt_tasks.front());
      next_task = this;
    }
    tbb::task_list at_list;
    for (auto && at : at_tasks) {
      using tbb_at_t = tbb_at_task<pc_type>; 
      auto * t = new(p->allocate_child()) tbb_at_t(std::move(at), _fc, _cih);
      at_list.push_back(*t);
    }
    if (dt_tasks.empty()) {  // we couldn't recycle anything
      set_ref_count(at_tasks.size() + 1);
      spawn_and_wait_for_all(at_list);
    } else {
      spawn(at_list);
    }
    return next_task;
  }
private:
  dt_type    _dt;
  fc_type &  _fc;
  cih_type  _cih;
};

}

#endif  // TBB_LL_HPP_

