#ifndef POISSON_H
#define POISSON_H

double poisson_log_partition(const double theta);

double poisson_log_partition_fderivative(const double theta);

void poisson_log_partition_derivatives(const double theta, double *b0,
                                       double *logb1, double *logb2);

#endif
