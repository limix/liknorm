#ifndef EXPFAM_H
#define EXPFAM_H

typedef
  void log_partition(double theta, double *b0, double *logb1, double *logb2);

typedef
  double log_partition0(double theta);

typedef
  double log_partition1(double theta);

enum lik_name;

typedef struct
{
  double         y;
  double         aphi;
  double         log_aphi;
  log_partition  *lp;
  log_partition0 *lp0;
  log_partition1 *lp1;
  double         lower_bound;
  double         upper_bound;
  enum lik_name  name;
} ExpFam;

#endif
