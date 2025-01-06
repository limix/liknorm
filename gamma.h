#ifndef PARTITION_GAMMA_H
#define PARTITION_GAMMA_H

double gamma_log_partition(const double theta);

double gamma_log_partition_fderivative(const double theta);

void gamma_log_partition_derivatives(const double theta, double *b0,
                                     double *logb1, double *logb2);

#endif
