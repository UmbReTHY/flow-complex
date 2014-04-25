#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include <stack>

#include "flow_complex.hpp"
#include "point_cloud.hpp"
#include "ascend_task.hpp"
#include "descend_task.hpp"
#include "utility.hpp"

namespace FC {

template <typename size_type,    // type that is capable of holding the indices
                                 // of the critical points
          bool Aligned = false,
          typename PointIterator,
          typename dim_type,
          typename num_threads_t,
          typename num_tol_t>
flow_complex<typename base_t<decltype((*PointIterator())[0])>::type,
             size_type>
compute_flow_complex (PointIterator begin, PointIterator end,
                      dim_type dim,
                      num_threads_t /*num_threads*/,
                      num_tol_t /*num_tolerance*/) {
  using number_type = typename base_t<decltype((*PointIterator())[0])>::type;
  using fc_type = flow_complex<number_type, size_type>;
  using pc_type = point_cloud<number_type, size_type, Aligned>;
  using at_type = ascend_task<pc_type>;
  using dt_type = descend_task<pc_type>;

  fc_type fc;
  // 1) init data structures
  pc_type pc(begin, end, dim);
  std::stack<at_type> qa;
  std::stack<dt_type> qd;
  // 2) seed initial ascend task(s)
  qa.emplace(pc);
  while (not qa.empty()) {
    at_type & at = qa.top();
    auto dh_handler = [&qd](dt_type && dt) {qd.push(std::move(dt));};
    // TODO define cp-Handler : takes cp and returns true if inserted, false otherwise
    // TODO define at-Handler: see dt Handler
    at.execute();
    qa.pop();
  }
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

