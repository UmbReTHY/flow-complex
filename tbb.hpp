#ifndef TBB_HPP_
#define TBB_HPP_

#include <vector>
#include <utility>
#include <memory>

#include <tbb/parallel_do.h>

#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"

namespace FC {

/**
  @brief Performs one step of an exploration of the maxima graph. In case a
         maximum is found, all descend tasks that can be executed from this
         maximum are also carried out.
*/
template <class point_cloud_type_>
class Task {
public:
  using pc_type = point_cloud_type_;
  using self_t = Task<pc_type>;
  using at_type = ascend_task<pc_type>;
  using dt_type = descend_task<pc_type>;
  using float_type = typename pc_type::number_type;
  using size_type = typename pc_type::size_type;
  using fc_type = flow_complex<float_type, size_type>;
  using item_t = std::shared_ptr<self_t>;
  using feeder_t = tbb::parallel_do_feeder<item_t>;
  
  static const int threshold = 0;

  Task(at_type at) {
    _atvec.push_back(std::move(at));
  }
  
  Task (Task const&) = delete;
  Task & operator=(Task const&) = delete;
  
  Task (Task &&) = default;
  Task & operator=(Task &&) = default;

  double weight() const {
    return _atvec.size() + _dtvec.size();
  }
  
  void add_task(at_type at) {
    _atvec.push_back(std::move(at));
  }
  
  void add_task(dt_type dt) {
    _dtvec.push_back(std::move(dt));
  }

  template <class ATCIHandler, class DTCIHandler>
  void execute(feeder_t & f, fc_type & fc,
               ATCIHandler & acih, DTCIHandler & dcih) {
    item_t share_task(new self_t());
    
    auto ath = [this, &share_task] (at_type at) {
      self_t * task = (weight() > threshold) ? share_task.get() : this;
      task->add_task(std::move(at));
    };
    auto dth = [this, &share_task] (dt_type dt) {
      self_t * task = (weight() > threshold) ? share_task.get() : this;
      task->add_task(std::move(dt));
    };

    work:
    while (true) {
      if (share_task->weight() > threshold) {
        f.add(share_task);
        share_task.reset(new self_t());
      }
      if (not _dtvec.empty())  // TODO measure depth-first / vs. breadt-first
        process_task(_dtvec, ath, dth, fc, dcih);
      else if (not _atvec.empty())
        process_task(_atvec, ath, dth, fc, acih);
      else
        break;
    }
    if (share_task->weight()) {
      std::swap(*share_task.get(), *this);
      goto work;
    }
  }
private:
  Task() = default;

  template <class task_t, class ... Params>
  void process_task(std::vector<task_t> & tvec, Params & ... params) {
    auto task(std::move(tvec.back()));
    tvec.pop_back();
    task.execute(std::forward<Params &>(params)...);
  }

  std::vector<at_type> _atvec;
  std::vector<dt_type> _dtvec;
};

}

#endif  // TBB_HPP_

