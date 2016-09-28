#ifndef DEFINITIONS_H
#define DEFINITIONS_H

typedef
  void log_partition(double theta, double *b0, double *logb1, double *logb2);

typedef
  double log_partition0(double theta);

typedef struct
{
  double         y;
  double         aphi;
  double         log_aphi;
  log_partition *lp;
  log_partition0 *lp0;
  double         left;
  double         right;
} ExpFam;

typedef struct
{
  double eta;
  double log_tau;
  double tau;
} Normal;

typedef struct
{
  double *log_zeroth; // log(mom0)
  double *u; // mom1/mom0
  double *v; // mom2/mom0
  int     n;
  double *A0;
  double *logA1;
  double *logA2;
  double *midiff;
} LikNormMachine;

#endif /* end of include guard: DEFINITIONS_H */
