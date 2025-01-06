#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

double liknorm_exponential_log_partition(const double theta);

double liknorm_exponential_log_partition_fderivative(const double theta);

void liknorm_exponential_log_partition_derivatives(const double theta,
                                                   double *b0, double *logb1,
                                                   double *logb2);

#endif
