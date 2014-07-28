#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <ostream>

namespace FC {

/**
  Logger to std::cout.
*/
class Logger {
public:
  typedef std::ostream& (*STRFUNC)(std::ostream&);

  template <typename Anything>
  Logger & operator<<(Anything const& at) {
    #ifdef LOGGING
    std::cout << at;
    #endif
    return *this;
  }
  
  /** to handle I/O manipulators, such as std::endl, correctly */
  Logger & operator<<(STRFUNC sf) {
    #ifdef LOGGING
    sf(std::cout);
    #endif
    return *this;
  }
};

}  // namespace FC

#endif  // LOGGER_HPP_

