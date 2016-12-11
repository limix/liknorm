#ifndef MACHINE_H
#define MACHINE_H

struct LikNormMachine
{
  double *log_zeroth; // log(mom0)
  double *u; // mom1/mom0
  double *v; // mom2/mom0
  int     size;
  double *A0;
  double *logA1;
  double *logA2;
  double *diff;
};

#endif
