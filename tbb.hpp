#ifndef TBB_HPP_
#define TBB_HPP_

#include <stack>
#include <utility>

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
class MaxTask {
public:
  using pc_type = point_cloud_type_;
  using at_type = ascend_task<pc_type>;

  MaxTask(at_type at) : _at(std::move(at)) {}

  template <class ATHandler,  class ATCIHandler, class DTCIHandler,
            typename float_t, typename size_type>
  void execute(ATHandler & ath, ATCIHandler & acih, DTCIHandler & dcih,
               flow_complex<float_t, size_type> & fc) {
    using dt_type = descend_task<pc_type>;
    std::stack<dt_type> dt_stack;
    auto dth = [&dt_stack] (dt_type dt) {dt_stack.push(std::move(dt));};
    _at.execute(ath, dth, fc, acih);
    while (not dt_stack.empty()) {
      auto dt(std::move(dt_stack.top()));
      dt_stack.pop();
      dt.execute(ath, dth, fc, dcih);
    }
  }
private:
  at_type _at;
};

}

#endif  // TBB_HPP_

