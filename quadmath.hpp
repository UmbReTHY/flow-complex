#ifndef QUADMATH_HPP
#define QUADMATH_HPP

#include <quadmath.h>
#include <ostream>

namespace std {
  __float128 abs(__float128 x) {
    return fabsq(x);
  }

  __float128 sqrt(__float128 x) {
    return sqrtq(x);
  }

}  // namespace std

#endif  // QUADMATH_HPP

