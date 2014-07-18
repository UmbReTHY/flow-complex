#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include <algorithm>
#include <utility>
#include <vector>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/task.h>

#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"
#include "logger.hpp"
#include "point_cloud.hpp"
#include "utility.hpp"

#include "tbb_ll.hpp"

namespace FC {

template <typename size_type,    // type that is capable of holding the indices
                                 // of the critical points
          bool Aligned = false,
          typename PointIterator,
          typename dim_type,
          typename num_threads_t>
flow_complex<typename base_t<decltype((*PointIterator())[0])>::type,
             size_type>
compute_flow_complex (PointIterator begin, PointIterator end,
                      dim_type dim,
                      num_threads_t /*num_threads*/ num_threads = 1) {
  Logger() << "*****************COMPUTE-START*************************\n";
  using number_type = typename base_t<decltype((*PointIterator())[0])>::type;
  using fc_type = flow_complex<number_type, size_type>;
  using pc_type = point_cloud<number_type, size_type, Aligned>;
  using at_type = ascend_task<pc_type>;
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
  
  // 3) seed initial ascend task(s)
  tbb::task_scheduler_init init(num_threads);
  tbb::task_list tl;
  for (num_threads_t t_id = 0; t_id < num_threads; ++t_id) {
    using tbb_at_t = tbb_at_task<pc_type>;
    auto * tbb_at = new(tbb::task::allocate_root()) tbb_at_t(at_type(pc), fc, acih);
    tl.push_back(*tbb_at);
  }
  // 4) process all tasks
  tbb::task::spawn_root_and_wait(tl);
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

