#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include <algorithm>
#include <stack>
#include <utility>
#include <vector>

#include <tbb/concurrent_unordered_set.h>

#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "flow_complex.hpp"
#include "logger.hpp"
#include "point_cloud.hpp"
#include "utility.hpp"

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
                      num_threads_t /*num_threads*/) {
  Logger() << "*****************COMPUTE-START*************************\n";
  using number_type = typename base_t<decltype((*PointIterator())[0])>::type;
  using fc_type = flow_complex<number_type, size_type>;
  using cp_type = typename fc_type::cp_type;
  using pc_type = point_cloud<number_type, size_type, Aligned>;
  using at_type = ascend_task<pc_type>;
  using dt_type = descend_task<pc_type>;
  using ci_type = circumsphere_ident<size_type>;

  // 1) init data structures
  pc_type pc(begin, end, dim);
  fc_type fc(dim, pc.size());
  std::stack<at_type>  qa;
  std::stack<dt_type>  qd;
  using ci_container = tbb::concurrent_unordered_set<ci_type, CIHash>;
  ci_container infproxy_cont;
  ci_container dci;
  // 2) create the handlers for task communication
  auto ath = [&qa] (at_type && at) {qa.push(std::move(at));};
  auto dth = [&qd] (dt_type && dt) {qd.push(std::move(dt));};
  auto cph = [&fc] (cp_type && cp) {return fc.insert(std::move(cp));};

  auto cih = [] (ci_container & ci_store, ci_type ci) {
    return ci_store.insert(ci).second;
  };
  auto acih = std::bind(cih, std::ref(infproxy_cont), std::placeholders::_1);
  auto dcih = std::bind(cih, std::ref(dci), std::placeholders::_1);
  
  // 3) seed initial ascend task(s)
  qa.emplace(pc);
  // 4) process all tasks - qa first = breadth first search
  while (not qa.empty() or not qd.empty()) {
    if (not qa.empty()) {
      at_type at(std::move(qa.top()));  // TODO code duplication: write "process task" function
      qa.pop();
      at.execute(dth, ath, fc, acih);
    } else {
      dt_type dt(std::move(qd.top()));
      qd.pop();  // TODO: these two steps above have to be atomic
      dt.execute(dth, ath, cph, dcih);
    }
  }
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

