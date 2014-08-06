#ifndef FLOW_COMPLEX_HPP_
#define FLOW_COMPLEX_HPP_

#include <cmath>

#include <algorithm>
#include <utility>
#include <map>
#include <ostream>
#include <istream>
#include <string>
#include <sstream>

#include <tbb/concurrent_unordered_set.h>

#include "critical_point.hpp"
#include "logger.hpp"

namespace FC {

template <typename nt, typename st> class flow_complex;
template <typename nt, typename st>
std::ostream & operator<<(std::ostream &, flow_complex<nt, st> const&);

template <typename _number_type, typename _size_type>
class flow_complex {
using self_t = flow_complex<_number_type, _size_type>;
  public:
    typedef _number_type                             number_type;
    typedef _size_type                                 size_type;
    typedef critical_point<_number_type, _size_type>     cp_type;

  private:
    using cp_container = tbb::concurrent_unordered_set<cp_type, CPHash>;
    
  public:
    typedef typename cp_container::iterator             iterator;
    typedef typename cp_container::const_iterator const_iterator;
  
    /**
      @brief initializes the flow complex with a maximum at infinity
    */
    flow_complex(size_type dim, size_type num_pts)
    : _minima(num_pts, nullptr) {
      _max_at_inf = &*_cps.insert(cp_type(dim)).first;
      // init fc with id-0 critical points
      for (size_type i = 0; i < num_pts; ++i)
        _minima[i] = insert(cp_type(&i, &i + 1, 0)).second;
    }
    // copy- and move-constructor
    flow_complex(flow_complex const&) = delete;
    flow_complex(flow_complex && tmp)
    : _max_at_inf(tmp._max_at_inf), _cps(), _minima(std::move(tmp._minima)) {
      _cps.swap(tmp._cps);
    }
    // copy- and move-assign
    flow_complex & operator=(flow_complex const&) = delete;
    flow_complex & operator=(flow_complex && rhs) {
      if (this != &rhs) {
        _max_at_inf = rhs._max_at_inf;
        _cps.swap(rhs._cps);
        _minima = std::move(rhs._minima);
      }
      return *this;
    }
  
    // iterators
    iterator begin() {
      return _cps.begin();
    }
    
    iterator end() {
      return _cps.end();
    }

    // const-iterators
    const_iterator begin() const {
      return _cps.begin();
    }

    const_iterator end() const {
      return _cps.end();
    }
    
    // modifiers
    void erase(cp_type const&);
    void erase(iterator);
    void erase(const_iterator);
    
    /**
      @return the first member of the pair is true, if the critical point was
              inserted by this call, and false if it already existed
    */
    std::pair<bool, cp_type *> insert(cp_type && cp) {
      std::pair<bool, cp_type *> r;
      auto ret_pair = _cps.insert(cp);
      r.first =   ret_pair.second;
      r.second = &*ret_pair.first;
      if (r.first)
        Logger() << "CRITICAL POINT IS NEW!\n";
      else
        Logger() << "CRITICAL POINT WAS FOUND ALREADY\n";
      return r;
    }
    
    cp_type const* max_at_inf() const {
      return _max_at_inf;
    }
    
    cp_type * minimum(size_type pos) const {
      return _minima[pos];
    }
    
    cp_type * find(cp_type const& cp) const {
      auto it = _cps.find(cp);
      return (it == _cps.end() ? nullptr : &*it);
    }
private:
  cp_type *              _max_at_inf;
  cp_container                  _cps;
  std::vector<cp_type *>     _minima;
    
  template <typename nt, typename st>
  friend std::istream & operator>>(std::istream &, flow_complex<nt, st> &);
};

template <typename number_type, typename size_type>
bool validate(flow_complex<number_type, size_type> const& fc) {
  std::map<size_type, int> hist;
  for (auto const& cp : fc)
    if (not cp.is_max_at_inf()) {
      auto r_pair = hist.emplace(cp.index(), 1);
      if (not r_pair.second)
        ++r_pair.first->second;
    }
  int sum = 0;
  for (auto const& el : hist)
    sum += (0 == (el.first % 2) ? el.second : -el.second);
  return 1 == sum;
}


template <typename nt, typename st>
std::ostream & operator<<(std::ostream & os, flow_complex<nt, st> const& fc) {
  for (auto const& cp : fc) {
    if (cp.is_max_at_inf()) {
      os << cp << cp.index() << std::endl;
    } else {
      os << cp;
      os << "| ";
      os << std::sqrt(cp.sq_dist());
      os << " ";
      for (auto it = cp.succ_begin(); it != cp.succ_end(); ++it)
        os << "| " << **it;
      os << std::endl;
    }
  }
  return os;
}

template <typename nt, typename st>
std::istream & operator>>(std::istream & is,
                          flow_complex<nt, st> & fc) {
  // types used
  using fc_type = flow_complex<nt, st>;
  using cp_type = typename fc_type::cp_type;
  using size_type = std::string::size_type;
  using float_t = typename fc_type::number_type;
  // init the given fc
  fc._cps.clear();
  fc._max_at_inf = nullptr;
  fc._minima.clear();
  // "global" variables for lambdas
  std::string line;
  std::vector<st> idxvec;
  float_t dist;
  auto const WS = " \t";  // WHITESPACE
  auto const DELIM = '|';
  // helper lambdas
  auto is_inf_line = [WS] (std::string const& s, size_type start = 0) {
    size_type pos1 = s.find_first_not_of(WS, start);
    size_type pos2 = s.find_first_of(WS, pos1);
    pos2 = (pos2 == std::string::npos) ? pos2 : pos2 - pos1;
    return s.substr(pos1, pos2) == "inf";
  };
  auto check_pos = [] (std::size_t pos) {
    if (pos == std::string::npos) throw std::invalid_argument("parse error");
  };
  auto parse_indices = [&] (std::string const& s, size_type start) {
    idxvec.clear();
    size_type end = s.find_first_of(DELIM, start);
    if (not is_inf_line(s, start)) {
      std::istringstream is(s.substr(start, end - start));
      st tmp;
      while (is >> tmp)
        idxvec.push_back(tmp);
    }
    return end;
  };
  auto parse_dist = [&dist] (std::string const& s, size_type start) {
    size_type end = s.find_first_of(DELIM, start);
    std::istringstream is(s.substr(start, end - start));
    is >> dist;
    return end;
  };
  // parsing loop
  while(std::getline(is, line)) {
    // skip empty and WS lines
    if (line.empty() or std::string::npos == line.find_first_not_of(WS))
      continue;
    if (is_inf_line(line)) {  // parse the cp at inf
      auto pos = line.find_first_not_of(WS);
      check_pos(pos);
      pos = line.find_first_of(WS, pos + 1);
      check_pos(pos);
      std::istringstream dim_stream(line.substr(pos));
      st dim;
      dim_stream >> dim;
      fc._max_at_inf = fc.insert(cp_type(dim)).second;
    } else {  // parse regular cps
      size_type start = 0;
      // parse indices of cp itself
      start = parse_indices(line, start);
      check_pos(start);  // there has to be a DELIM
      // parse distance
      start = parse_dist(line, ++start);
      auto rpair = fc.insert(cp_type(idxvec.begin(), idxvec.end(), 0));
      auto * cp_ptr = rpair.second;
      // update _minima array
      if (0 == cp_ptr->index()) {
        st const min_pos = *cp_ptr->idx_begin();
        if (fc._minima.size() <= min_pos)
          fc._minima.resize(min_pos + 1);
        fc._minima[min_pos] = cp_ptr;
      }
      cp_ptr->_sq_dist = dist * dist;  // even if the insert was new
      // parse succsessors
      while (std::string::npos != start) {
        start = parse_indices(line, ++start);
        if (idxvec.empty()) {  // succ is cp at inf
          // assumes the data format always lists the cp at inf first
          assert(fc._max_at_inf);
          cp_ptr->add_successor(fc._max_at_inf);
        } else {
          auto rpair = fc.insert(cp_type(idxvec.begin(), idxvec.end(), 0));
          auto const* succ_ptr = rpair.second;
          cp_ptr->add_successor(succ_ptr);
        }
      }
    }
  }
  return is;
}

}  // namespace FC

#endif  // FLOW_COMPLEX_HPP_

