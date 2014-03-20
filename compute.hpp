#ifndef COMPUTE_HPP_
#define COMPUTE_HPP_

#include "flow_complex.hpp"
#include "flow_point.hpp"
#include "point_cloud.hpp"

namespace FC {

template <typename number_type, typename size_type, typename num_threads_t>
flow_complex<number_type, size_type>
compute_flow_complex (number_type const** begin, number_type const** end,
                      size_type dim,
                      num_threads_t /*num_threads*/,
                      number_type /*num_tolerance*/) {  
  flow_complex<number_type, size_type> fc;
  
  // 0) form the point cloud object
  point_cloud<number_type, size_type> pc(begin, end, dim);
  // 1) find first maximum
  flow_point<number_type, size_type> fp(pc);
  while (not fp.is_proxy_at_inf() and not fp.is_finite_max())
    fp.ascend();
  
  return fc;
}

}  // namespace FC

#endif  // COMPUTE_HPP_

