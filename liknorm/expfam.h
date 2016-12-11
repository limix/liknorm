#ifndef EXPFAM_H
#define EXPFAM_H

typedef void log_partition_derivatives(const double theta, double *b0,
                                       double *logb1, double *logb2);

typedef double log_partition_fderivative(const double theta);

enum lik_name;

typedef struct {
  double y;
  double aphi;
  double log_aphi;
  log_partition_derivatives *lpd;
  log_partition_fderivative *lpfd;
  double lower_bound;
  double upper_bound;
  enum lik_name name;
} ExpFam;

#endif
