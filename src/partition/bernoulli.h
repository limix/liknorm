#ifndef BERNOULLI_H
#define BERNOULLI_H

#include "partition/binomial.h"

static inline double bernoulli_log_partition(const double theta)
{
    return binomial_log_partition(theta);
}

static inline double bernoulli_log_partition_fderivative(const double theta)
{
    return binomial_log_partition_fderivative(theta);
}

static inline void bernoulli_log_partition_derivatives(const double theta, double *b0,
                                                       double *logb1, double *logb2)
{
    binomial_log_partition_derivatives(theta, b0, logb1, logb2);
}

#endif
