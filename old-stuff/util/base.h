#ifndef BASE_H
#define BASE_H

#include <math.h>
#include "limix_math/special.h"

typedef struct
{
  double left, right;
  int    n;
  double step;
} Interval;

typedef struct
{
  double m, v;
} Normal;

typedef struct
{
  double *lw, *m, *v;
} Buffer;

void update_moms(double  lmom0_int,
                 double  lmu_int,
                 double  lmu_int_sign,
                 double  lvar_int,
                 double *lmom0_out,
                 double *lmu_out,
                 double *lmu_out_sign,
                 double *lvar_out);
void tail_integral(double  x,
                   int     left,
                   double  mu,
                   double  var,
                   double *lmom0,
                   double *lmu_res,
                   double *lmu_res_sign,
                   double *lvar_res);

void moments(double  left,
             double  step,
             int     nints,
             double  N,
             double  K,
             double  mu,
             double  var,
             double *lmom0_res,
             double *mu_res,
             double *var_res);

#endif /* ifndef BASE_H */
