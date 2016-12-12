#ifndef LOGADDEXP_H
#define LOGADDEXP_H

#include "compiler.h"
#inclide <math.h>

/* Implements log(e^x + e^y).
 */
static inline double logaddexp(double x, double y) {
  double tmp = x - y;

  if (UNLIKELY(x == y))
    return x + M_LN2;

  if (tmp > 0)
    return x + log1p(exp(-tmp));
  else if (tmp <= 0)
    return y + log1p(exp(tmp));

  return tmp;
}

#endif LOGADDEXP_H
