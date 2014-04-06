#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include "flow_complex.hpp"
#include "point_cloud.hpp"
#include "ascend_task.hpp"

namespace FC {

// TODO enable template flag that enables the user to specify
//      if the data pointers are SSE aligned, or not
template <typename number_type, typename size_type, typename num_threads_t>
flow_complex<number_type, size_type>
// TODO generalize so that any iterators can be passed
compute_flow_complex (number_type const** begin, number_type const** end,
                      size_type dim,
                      num_threads_t /*num_threads*/,
                      number_type /*num_tolerance*/) {
  using fc_type = flow_complex<number_type, size_type>;
  using pc_type = point_cloud<number_type, size_type>;
  using at_type = ascend_task<number_type, size_type>;
  fc_type fc;
  
  // 0) form the point cloud object
  pc_type pc(begin, end, dim);
  // 1) find first maximum
  at_type at(pc);
  at.execute();
  // TODO enlist maximum (inf or regular)
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

