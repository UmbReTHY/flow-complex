#ifndef CLEAN_HPP_
#define CLEAN_HPP_

#include <algorithm>
#include <functional>
#include <iterator>
#include <map>
#include <vector>

#include "flow_complex.hpp"

namespace FC {

template <typename nt, typename st>
bool exists_path(critical_point<nt, st> const* from,
                 critical_point<nt, st> const* to) {
  using namespace std::placeholders;
  auto exists_path_to = std::bind(&exists_path<nt, st>, _1, to);
  return to->index() > from->index() and
         ((from->succ_end() !=
          std::find(from->succ_begin(), from->succ_end(), to)) or
          (from->succ_end() !=
          std::find_if(from->succ_begin(), from->succ_end(), exists_path_to))
         );
}

/**
@brief Let x, y be predecessors of a critical point z. Because of the way we
       explore the Hasse Diagram, it may happen that x is also the predecessor
       of y, or the predecessor of one of y's predecessors. In that case, we
       remove the incidence from x to z. This function does
*/
template <typename nt, typename st>
flow_complex<nt, st> clean_incidences(flow_complex<nt, st> fc) {
  using fc_type = flow_complex<nt, st>;
  using cp_type = typename fc_type::cp_type;

  // contains a pair for every cp, that has redundant successors
  std::map<cp_type *, std::vector<cp_type const*>> to_delete;
  // detect redundant successors
  for (auto & cp : fc) {  // TODO parallelize this loop
    std::vector<cp_type const*> redundant_succs;
    for (auto from_it = cp.succ_begin(); from_it != cp.succ_end(); ++from_it)
      for (auto to_it = cp.succ_begin(); to_it != cp.succ_end(); ++to_it)
        if (exists_path(*from_it, *to_it))
          redundant_succs.push_back(*to_it);
    if (not redundant_succs.empty())
      to_delete.emplace(&cp, std::move(redundant_succs));
  }
  // clean the flow complex from redundant incidences
  for (auto const& p : to_delete)
    for (auto succ : p.second)
      p.first->erase(succ);
  return fc;
}

}  // namespace FC

#endif  // CLEAN_HPP_

