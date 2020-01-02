#include "partition/gamma.h"
#include "hide.h"
#include <float.h>
#include <math.h>

HIDE double gamma_log_partition(const double theta)
{
    return -log(fmax(DBL_EPSILON, -theta));
}

HIDE double gamma_log_partition_fderivative(const double theta)
{
    return -log(fmax(DBL_EPSILON, -theta));
}

HIDE void gamma_log_partition_derivatives(const double theta, double *b0, double *logb1,
                                          double *logb2)
{
    *b0 = -log(fmax(DBL_EPSILON, -theta));
    *logb1 = *b0;
    *logb2 = 2 * (*b0);
}
