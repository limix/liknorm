#include "partition/exponential.h"
#include "hide.h"
#include <assert.h>
#include <float.h>
#include <math.h>

HIDE double exponential_log_partition(const double theta)
{
    return -log(fmax(DBL_EPSILON, -theta));
}

HIDE double exponential_log_partition_fderivative(const double theta)
{
    return -log(fmax(DBL_EPSILON, -theta));
}

HIDE void exponential_log_partition_derivatives(const double theta, double *b0,
                                                double *logb1, double *logb2)
{
    *b0 = -log(fmax(DBL_EPSILON, -theta));
    *logb1 = *b0;
    *logb2 = 2 * (*b0);
}
