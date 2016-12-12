#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

double exponential_log_partition(const double theta);

double exponential_log_partition_fderivative(const double theta);

void exponential_log_partition_derivatives(const double theta, double *b0,
                                           double *logb1, double *logb2);

#endif
