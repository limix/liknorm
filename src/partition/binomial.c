#include "partition/binomial.h"
#include <float.h>
#include <math.h>

double binomial_log_partition(const double theta)
{
    return theta < log(DBL_MAX) ? log1p(exp(theta)) : theta;
}

double binomial_log_partition_fderivative(const double theta)
{
    return theta < log(DBL_MAX) ? theta - log1p(exp(theta)) : 0;
}

void binomial_log_partition_derivatives(const double theta, double *b0, double *logb1,
                                        double *logb2)
{
    if (theta < log(DBL_MAX)) {
        double log1p_ = log1p(exp(theta));
        *b0 = log1p_;
        *logb1 = theta - log1p_;
        *logb2 = theta - 2 * log1p_;
    } else {
        *b0 = theta;
        *logb1 = 0;
        *logb2 = -theta;
    }
}
