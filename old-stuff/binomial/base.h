#ifndef BASE_H
#define BASE_H

#include <math.h>
#include "limix_math/special.h"

double log_ulike_prior(double x, double N, double K, double mu, double var);
void update_moms(double lmom0_int, double lmu_int, double lmu_int_sign,
                 double lvar_int, double* lmom0_out, double* lmu_out,
                 double* lmu_out_sign, double* lvar_out);
void tail_integral(double x, int left, double mu, double var,
                  double* lmom0, double* lmu_res, double* lmu_res_sign,
                  double* lvar_res);
void get_extrema(double N, double K, double mu, double var,
                 double* left, double* right);
static inline double get_mode(double N, double K)
{
  if (N == K)
    return lmath_normal_icdf(1./pow(2, 1/N));
  if (K == 0)
    return -lmath_normal_icdf(1./pow(2, 1/N));
  return lmath_normal_icdf(K / N);
}
void integrate_window(double l_, double r_,
                      double N, double K,
                      double mu, double var,
                      double* lmom0, double* lmu_res, double* lmu_res_sign,
                      double* lvar_res);
void moments(double left, double step, int nints,
             double N, double K, double mu, double var,
             double* lmom0_res, double* mu_res, double* var_res);

#endif
