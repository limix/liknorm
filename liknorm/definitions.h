#ifndef DEFINITIONS_H
#define DEFINITIONS_H

typedef struct
{
  double *log_zeroth; // log(mom0)
  double *u; // mom1/mom0
  double *v; // mom2/mom0
  int     n;
  double *A0;
  double *logA1;
  double *logA2;
  double *diff;
} LikNormMachine;

#endif /* end of include guard: DEFINITIONS_H */
