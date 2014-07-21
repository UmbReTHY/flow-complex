#ifndef MAKROS_H_
#define MAKROS_H_

#if __GNUC__ >= 3
  # define __pure         __attribute__ ((pure))
  # define __const        __attribute__ ((const))
  # define likely(x)      __builtin_expect (!!(x), 1)
  # define unlikely(x)    __builtin_expect (!!(x), 0)
#else
  # define _  _pure         /* no pure */
  # define _  _const        /* no const */
  # define likely(x)      (x)
  # define unlikely(x)    (x)
#endif

#endif  // MAKROS_H_

