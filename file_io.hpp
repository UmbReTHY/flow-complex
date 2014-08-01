#ifndef FILE_IO_HPP_
#define FILE_IO_HPP_

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <ostream>

namespace FC {

template <typename float_t>
class point_store {
public:
  typedef std::vector<float_t>                          point_t;
  typedef std::size_t                                   size_type;
  typedef typename std::vector<point_t>::iterator       iterator;
  typedef typename std::vector<point_t>::const_iterator const_iterator;

  point_store() = default;
  point_store(point_store const&) = default;
  point_store(point_store &&    ) = default;

  point_store & operator=(point_store const&) = default;
  point_store & operator=(point_store &&    ) = default;
  
  iterator begin() {return _pts.begin();}
  iterator end() {return _pts.end();}
  const_iterator cbegin() const {return _pts.cbegin();}
  const_iterator cend() const {return _pts.cend();}
  
  point_t & operator[](size_type pos) {return _pts[pos];}
  point_t const& operator[](size_type pos) const {return _pts[pos];}
  
  size_type size() const {return _pts.size();}
  
  point_t & add_point(size_type dim) {_pts.emplace_back(dim); return _pts.back();}
private:
  std::vector<point_t> _pts;
  size_type _dim;
};

template <typename float_t>
std::ostream & operator<<(std::ostream & os, point_store<float_t> const& ps) {
  for (auto it = ps.cbegin(); it != ps.cend(); ++it) {
    auto & p = *it;
    for (auto e : p)
      os << e << " ";
    os << "\n";
  }
  return os;
}

template <typename float_t>
point_store<float_t> from_file(char const* filename) {
  std::ifstream f(filename);
  if (not f)
    throw std::runtime_error("could not open file");
  std::string line;
  point_store<float_t> pts;
  std::vector<float_t> tmp_p;
  while (std::getline(f, line)) {
    if (line.empty())
      continue;
    std::istringstream linestream(line);
    tmp_p.clear();
    float_t tmpf;
    while (linestream >> tmpf)
      tmp_p.push_back(tmpf);
    if (tmp_p.size()) {
      auto p_it = pts.add_point(tmp_p.size()).begin();
      for (auto e : tmp_p)
        *p_it++ = e;
    }
  }
  return pts;
}

}

#endif  // FILE_IO_HPP_

