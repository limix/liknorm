#include "partition/poisson.h"
#include <float.h>
#include <math.h>

double poisson_log_partition(const double theta)
{
    return exp(fmin(theta, log(DBL_MAX)));
}

double poisson_log_partition_fderivative(const double theta) { return theta; }

void poisson_log_partition_derivatives(const double theta, double *b0, double *logb1,
                                       double *logb2)
{
    *b0 = exp(fmin(theta, log(DBL_MAX)));
    *logb2 = *logb1 = theta;
}
