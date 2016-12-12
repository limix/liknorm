#ifndef BINOMIAL_H
#define BINOMIAL_H

double binomial_log_partition(const double theta);

double binomial_log_partition_fderivative(const double theta);

void binomial_log_partition_derivatives(const double theta, double *b0,
                                        double *logb1, double *logb2);

#endif
