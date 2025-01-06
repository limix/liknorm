#include "geometric.h"
#include <float.h>
#include <math.h>

double liknorm_geometric_log_partition(const double theta)
{
    return -log1p(-exp(theta));
}

double liknorm_geometric_log_partition_fderivative(const double theta)
{
    return theta - log1p(-exp(theta));
}

void liknorm_geometric_log_partition_derivatives(const double theta, double *b0,
                                         double *logb1, double *logb2)
{
    const double log1p_ = -log1p(-exp(theta));
    *b0 = log1p_;
    *logb1 = theta + log1p_;
    *logb2 = theta + 2 * log1p_;
}
