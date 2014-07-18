#ifndef TBB_HPP_
#define TBB_HPP_

#include <functional>
#include <utility>

#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"

namespace FC {

template <class point_cloud_t>
class abstract_task {
public:
  typedef point_cloud_t                            pc_type;
  typedef typename pc_type::number_type        number_type;
  typedef typename pc_type::size_type            size_type;
  typedef flow_complex<number_type, size_type>     fc_type;
  typedef ascend_task<pc_type>                     at_type;
  typedef descend_task<pc_type>                    dt_type;
  typedef std::function<void(at_type &&)>         ath_type;
  typedef std::function<void(dt_type &&)>         dth_type;

  virtual void execute(ath_type ath, dth_type dth, fc_type & fc) = 0;
  
  virtual ~abstract_task() {
  }
};

template <class task_t, class cih_type>
class tbb_task_wrapper : public abstract_task<typename task_t::point_cloud_type> {
public:
  typedef typename task_t::point_cloud_type  pc_type;
  typedef abstract_task<pc_type>              base_t;
  typedef typename base_t::ath_type         ath_type;
  typedef typename base_t::dth_type         dth_type;
  typedef typename base_t::fc_type           fc_type;

  tbb_task_wrapper(task_t task,  cih_type & cih)
  : _task(std::move(task)), _cih(cih) {
  }

  virtual void execute(ath_type ath, dth_type dth, fc_type & fc) {
    _task.execute(ath, dth, fc, _cih);
  }

private:
  task_t     _task;
  cih_type &  _cih;
};

template <class task_t, class cih_t>
abstract_task<typename task_t::point_cloud_type> *
make_task_wrapper(task_t task, cih_t & cih) {
  using wrapper_t = tbb_task_wrapper<task_t, cih_t>;
  return new wrapper_t(std::move(task), cih);
}

}

#endif  // TBB_HPP_

