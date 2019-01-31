#ifndef NBINOMIAL_H
#define NBINOMIAL_H

double nbinomial_log_partition(const double theta);

double nbinomial_log_partition_fderivative(const double theta);

void nbinomial_log_partition_derivatives(const double theta, double *b0,
                                         double *logb1, double *logb2);

#endif
