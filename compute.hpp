#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include "flow_complex.hpp"
#include "point_cloud.hpp"
#include "ascend_task.hpp"
#include "utility.hpp"

namespace FC {

template <typename size_type,    // type that is capable of holding the indices
                                 // of the critical points
          bool Aligned = false,
          typename PointIterator,
          typename dim_type,
          typename num_tol_t>
flow_complex<typename base_t<decltype(**PointIterator())>::type,
             size_type>
compute_flow_complex (PointIterator begin, PointIterator end,
                      dim_type dim,
                      int /*num_threads*/,
                      num_tol_t /*num_tolerance*/) {
  using number_type = typename base_t<decltype(**PointIterator())>::type;
  using fc_type = flow_complex<number_type, size_type>;
  using pc_type = point_cloud<number_type, size_type, Aligned>;
  using at_type = ascend_task<pc_type>;

  fc_type fc;
  // 0) form the point cloud object
  pc_type pc(begin, end, dim);
  // 1) seed initial ascend task(s)
  at_type at(pc);
  // TODO function/object coordinating construction
  //      manages task queues and the flow complex structure
  //      tasks get a reference to this object on construction
  at.execute();
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

