#ifndef LOGADDEXP_H
#define LOGADDEXP_H

#include "compiler.h"
#include <math.h>

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

static inline double logaddexp_array(double *x, int n, double xmax) {
  double total = 0;

  for (int i = 0; i < n; ++i)
    total += exp(x[i] - xmax);

  return xmax + log(total);
}

#endif
