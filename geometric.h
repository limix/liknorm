#ifndef LIKNORM_GEOMETRIC_H
#define LIKNORM_GEOMETRIC_H

double liknorm_geometric_log_partition(const double theta);

double liknorm_geometric_log_partition_fderivative(const double theta);

void liknorm_geometric_log_partition_derivatives(const double theta, double *b0,
                                                 double *logb1, double *logb2);

#endif
