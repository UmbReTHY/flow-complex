#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include <algorithm>
#include <utility>
#include <vector>

#include <glog/logging.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_do.h>

#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"
#include "point_cloud.hpp"
#include "utility.hpp"
#include "tbb.hpp"

namespace FC {

/**
  @num_threads if negative, use Intel TBB default
*/
template <typename size_type,    // type that is capable of holding the indices
                                 // of the critical points
          bool Aligned = false,
          typename PointIterator,
          typename dim_type>
flow_complex<typename base_t<decltype((*PointIterator())[0])>::type,
             size_type>
compute_flow_complex (PointIterator begin, PointIterator end,
                      dim_type dim, int num_threads) {
  DLOG(INFO) << "*****************COMPUTE-START*************************\n";
  using number_type = typename base_t<decltype((*PointIterator())[0])>::type;
  using fc_type = flow_complex<number_type, size_type>;
  using pc_type = point_cloud<number_type, size_type, Aligned>;
  using at_type = ascend_task<pc_type>;
  using dt_type = descend_task<pc_type>;
  using ci_type = circumsphere_ident<size_type>;

  // 1) init data structures
  pc_type pc(begin, end, dim);
  fc_type fc(dim, pc.size());

  using ci_container = tbb::concurrent_unordered_set<ci_type, CIHash>;
  ci_container infproxy_cont;
  ci_container dci;
  // 2) create the handlers for task communication
  auto cih = [] (ci_container & ci_store, ci_type ci) {
    return ci_store.insert(ci).second;
  };
  auto acih = std::bind(cih, std::ref(infproxy_cont), std::placeholders::_1);
  auto dcih = std::bind(cih, std::ref(dci), std::placeholders::_1);
  
  // 3) seed initial ascend task(s)
  using item_t = std::shared_ptr<abstract_task<pc_type>>;
  using feeder_t = tbb::parallel_do_feeder<item_t>;
  if (num_threads < 0) num_threads = tbb::task_scheduler_init::default_num_threads();
  std::vector<item_t> tasks;
  tbb::task_scheduler_init init(num_threads);
  for (int i = 1; i <= num_threads; ++i)
    tasks.push_back(item_t(make_task_wrapper(at_type(pc), acih)));
  // 4) process all tasks
  tbb::parallel_do(tasks.begin(), tasks.end(), [&](item_t item, feeder_t & f) {
    using ttb_at_t = tbb_task_wrapper<at_type, decltype(acih)>;
    using ttb_dt_t = tbb_task_wrapper<dt_type, decltype(dcih)>;
    auto ath = [&] (at_type && at) {f.add(item_t(new ttb_at_t(std::move(at), acih)));};
    auto dth = [&] (dt_type && dt) {f.add(item_t(new ttb_dt_t(std::move(dt), dcih)));};
    item->execute(ath, dth, fc);
  });
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

