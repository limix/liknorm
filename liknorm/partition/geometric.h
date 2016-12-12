#ifndef GEOMETRIC_H
#define GEOMETRIC_H

double geometric_log_partition(const double theta);

double geometric_log_partition_fderivative(const double theta);

void geometric_log_partition_derivatives(const double theta, double *b0,
                                         double *logb1, double *logb2);

#endif
