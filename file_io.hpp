#ifndef FILE_IO_HPP_
#define FILE_IO_HPP_

#include <fstream>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <glog/logging.h>

#include "utility.hpp"

namespace FC {

template <typename number_type_, typename size_type_>
class point_store {
public:
  typedef number_type_                                     number_type;
  typedef std::vector<number_type_>                        point_type;
  typedef size_type_                                       size_type;
  typedef typename std::vector<point_type>::iterator       iterator;
  typedef typename std::vector<point_type>::const_iterator const_iterator;

  point_store() = default;
  point_store(point_store const&) = default;
  point_store(point_store &&    ) = default;

  point_store & operator=(point_store const&) = default;
  point_store & operator=(point_store &&    ) = default;
  
  iterator begin() {return pts_.begin();}
  iterator end() {return pts_.end();}
  const_iterator cbegin() const {return pts_.cbegin();}
  const_iterator cend() const {return pts_.cend();}
  
  point_type& operator[](size_type pos) {return pts_[pos];}
  const point_type& operator[](size_type pos) const {return pts_[pos];}
  
  size_type dim() const {
    CHECK(size() > 0);
    return convertSafelyTo<size_type>(pts_[0].size());
  }

  size_type size() const { return convertSafelyTo<size_type>(pts_.size());}
  
  // not thread-safe
  point_type& add_point(size_type dim) {
    if (size() > 0) {CHECK(dim == this->dim());}
    pts_.emplace_back(dim);
    return pts_.back();
  }
private:
  std::vector<point_type> pts_;
};

template <typename number_type, typename size_type>
std::ostream & operator<<(std::ostream & os,
                          const point_store<number_type, size_type> & ps) {
  for (auto it = ps.cbegin(); it != ps.cend(); ++it) {
    auto & p = *it;
    for (auto e : p)
      os << std::setprecision(std::numeric_limits<float_t>::digits10)
         << e << " ";
    os << "\n";
  }
  return os;
}

template <typename number_type, typename size_type>
std::istream & operator>>(std::istream & is,
                          point_store<number_type, size_type> & pts) {
  std::string line;
  std::vector<number_type> tmp_p;
  while (std::getline(is, line)) {
    if (line.empty())
      continue;
    std::istringstream linestream(line);
    tmp_p.clear();
    double tmpf;
    while (linestream >> tmpf)
      tmp_p.push_back(tmpf);
    if (tmp_p.size()) {
      const size_type dim = convertSafelyTo<size_type>(tmp_p.size());
      auto p_it = pts.add_point(dim).begin();
      for (auto e : tmp_p)
        *p_it++ = e;
    }
  }
  return is;
}

}  // namespace FC

#endif  // FILE_IO_HPP_

