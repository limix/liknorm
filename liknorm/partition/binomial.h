#ifndef BINOMIAL_H
#define BINOMIAL_H

#include <float.h>
#include <math.h>

extern const double max_exp;

static inline double binomial_log_partition_fderivative(const double theta)
{
  return -theta < max_exp ? -log1p(exp(-theta)) : theta;
}

static inline void binomial_log_partition_derivatives(
  const double theta, double *b0, double *logb1, double *logb2)
{
  double log1p_;
  if (-theta < max_exp)
  {
    log1p_ = log1p(exp(-theta));
    *b0 = theta + log1p_;
    *logb1 = -log1p_;
    *logb2 = - theta - 2 * log1p_;
  }
  else
  {
    *b0 = 0;
    *logb1 = theta;
    *logb2 = - theta - 2 * log1p_;
  }
}

#endif
