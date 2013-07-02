#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <tbb/concurrent_unordered_set.h>

#include "critical_point.hpp"

namespace FC {

template <typename number_type>
class flow_complex {
  using cp_container = tbb::concurrent_unordered_set<critical_point>;

  public:
    flow_complex()



  private:
     _critical_points;
};

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_
