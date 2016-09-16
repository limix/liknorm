#include <stdio.h>
#include <math.h>
#include "base.h"

/*
 # Definitions
   - Exponential-family distribution (likelihood):
   f(y|\eta) = h(y) \exp(\eta \cdot T(y) - A(\eta))

   - Log-partition function: A(\eta))

   - ulikelihood:
   l_u(y|\eta) = \exp(\eta \cdot T(y) - A(\eta))

 # User input
   - Values h(y) and T(y).
   - Functions A(\cdot), A'(\cdot), and A''(\cdot).
 */
inline
double logsumexp(double *x, int n)
{
  double r = x[0];
  int    i;

  for (i = 1; i < n; i++) r = lmath_logaddexp(r, x[i]);
  return r;
}

inline double normal_logpdf2(double x, double mu, double var)
{
  return lmath_normal_logpdf((x - mu) / sqrt(var)) - log(var) / 2.0;
}

// Computes \log l_u(y|\eta) + \log \mathcal{N}(\eta | m, v)
inline
double log_ulike_prior(double Ty, double (*A)(double, void *), void *args,
                       Normal *normal, double eta)
{
  return eta * Ty - (*A)(eta, args) + normal_logpdf2(eta, normal->m,
                                                     normal->v);
}

// Shrink the interval based on the upper bound on the AUC estimation.
int shrink_interval(double (*f)(double, void *), void *args,
                    Interval *interval)
{
  // log_limit is actually around two times it but
  // the reason I use half is that I know that at
  // the end I will do lw[i] - logsumexp(lw).
  double log_limit = -344;
  Interval *i      = interval;

  while ((*f)(i->left, args) <= log_limit) i->left += i->step;

  while ((*f)(i->right, args) <= log_limit) i->right -= i->step;

  // error report
  if (i->left >= i->right) return 1;

  i->step = (i->right - i->left) / i->n;

  return 0;
}

// // Compute log_ulike_prior for each point in the interval boundaries.
// void log_curve_ordinates(double (*lulik)(double x), void *args,
//                          Interval *interval, Normal *normal,
//                          double *log_ordinates)
// {
//   int i;
//   double x = interval->left;
//
//   for (i = 0; i < interval->n + 1; i++)
//   {
//     log_ordinates[i] = log_ulike_prior(lulik, args, normal, x);
//     x               += interval->step;
//   }
// }

// void shrink_intervals(double (*lulik)(double *x),
//                       void *args,
//                       Interval *interval,
//                       Normal *normal,
//                       double *heights, double *w)
//
//
// // Having the heights, it remains to area under the curve of each interval.
// double compute_weights(Interval * interval, Normal * normal,
//                        double * heights, double * w)
// {
//   double hl, hr, c;
//   double r, llc;
//   double acum;
//   int    i;
//
//   for (i = 0; i < interval->n; i++)
//   {
//     r = interval->left + interval->step;
//
//     hl = heights[i];
//     hr = heights[i + 1];
//     c  = (hr - hl) / interval->step;
//
//     llc = log(fabs(c));
//
//     if (c > 0)
//     {
//       w[i] = hl - interval->left * c + lmath_logaddexps(r * c,
//                                                         interval->left * c,
//                                                         1,
//                                                         -1) - llc;
//     } else
//     {
//       w[i] = hl - interval->left * c + lmath_logaddexps(r * c,
//                                                         interval->left * c,
//                                                         -1,
//                                                         1) - llc;
//     }
//
//     if (i == 0) acum = w[i];
//     else acum = lmath_logaddexp(acum, w[i]);
//
//     interval->left = r;
//   }
//
//   return acum;
// }

// void update_moms(double lmom0_int, double lmu_int, double lmu_int_sign,
//                  double lvar_int, double *lmom0_out, double *lmu_out,
//                  double *lmu_out_sign, double *lvar_out)
// {
//   double sign;
//
//   *lmom0_out = lmath_logaddexp(*lmom0_out, lmom0_int);
//
//   *lmu_out = lmath_logaddexpss(*lmu_out, lmu_int + lmom0_int, *lmu_out_sign,
//                                lmu_int_sign, &sign);
//   *lmu_out_sign = sign;
//
//   *lvar_out = lmath_logaddexp(*lvar_out, lvar_int + lmom0_int);
//   *lvar_out = lmath_logaddexp(*lvar_out, 2 * lmu_int + lmom0_int);
// }

// void tail_integral(double x, int left, double mu, double var,
//                    double *lmom0, double *lmu_res, double *lmu_res_sign,
//                    double *lvar_res)
// {
//   double sig  = sqrt(var);
//   double lsig = log(sig);
//   double m    = (x - mu) / sig;
//   double lpdf = lmath_normal_logpdf(m);
//
//   if (left) *lmom0 = lmath_normal_logsf(m);
//   else *lmom0 = lmath_normal_logcdf(m);
//
//   double sign;
//
//   if (left) sign = +1.0;
//   else sign = -1.0;
//
//   *lmu_res = mu + sign * exp(lpdf - *lmom0 + lsig);
//
//   if (*lmu_res >= 0.) *lmu_res_sign = +1.0;
//   else
//   {
//     *lmu_res_sign = -1.0;
//     *lmu_res     *= -1.0;
//   }
//   *lmu_res = log(*lmu_res);
//
//   *lvar_res = var + sign * sig * (mu + x) * exp(lpdf - *lmom0)
//               - exp(2 * *lmu_res) + mu * mu;
//   *lvar_res = log(*lvar_res);
// }

void integrate_window(void (*f_fl_fll)(double,
                                       void *,
                                       double *,
                                       double *,
                                       double *), void *args,
                      double l_, double r_,
                      Normal *normal,
                      double *lw, double *m, double *v)
{
  double x0 = l_ + (r_ - l_) / 2;
  double a, b, c;

  (*f_fl_fll)(x0, args, &a, &b, &c);

  double mu    = normal->m;
  double var   = normal->v;
  double alpha = a - b * x0 + c * x0 * x0 / 2 - mu * mu / (2 * var);
  double beta  = b - c * x0 + mu / var;
  double gamma = c / 2 - 1 / (2 * var);

  double mu_  = -beta / (2 * gamma);
  double var_ = -1 / (2 * gamma);
  double sig_ = sqrt(var_);

  double lcoef = -beta * beta / (4 * gamma) + alpha;
  lcoef -= log(-2 * gamma * var) / 2;

  double ml = (l_ - mu_) / sig_;
  double mr = (r_ - mu_) / sig_;

  double lpdf_r = lmath_normal_logpdf(mr);
  double lpdf_l = lmath_normal_logpdf(ml);
  double lcdf_l, lcdf_r;
  double lcdf_diff, lpdf_diff;
  double lsf_l, lsf_r;

  if (ml + mr >= 0)
  {
    lsf_l     = lmath_normal_logsf(ml);
    lsf_r     = lmath_normal_logsf(mr);
    lcdf_diff = lmath_logaddexps(lsf_l, lsf_r, 1, -1);
  }
  else
  {
    lcdf_l    = lmath_normal_logcdf(ml);
    lcdf_r    = lmath_normal_logcdf(mr);
    lcdf_diff = lmath_logaddexps(lcdf_r, lcdf_l, 1, -1);
  }

  *lmom0 = lcoef + lcdf_diff;
  double lpdf_diff_sign;

  if (lpdf_l == lpdf_r)
  {
    *lmu_res = log(fabs(mu_));

    if (mu_ >= 0) *lmu_res_sign = +1;
    else *lmu_res_sign = -1;
  }
  else
  {
    lpdf_diff = lmath_logaddexpss(lpdf_l, lpdf_r, +1, -1, &lpdf_diff_sign);
    *lmu_res  = mu_ + sig_ * exp(lpdf_diff - lcdf_diff) * lpdf_diff_sign;

    if (*lmu_res >= 0) *lmu_res_sign = +1;
    else *lmu_res_sign = -1;
    *lmu_res = log(fabs(*lmu_res));
  }
  double t1_sign;
  double t1 = lmath_logaddexpss(lpdf_l, lpdf_r, ml, -mr, &t1_sign);
  t1 = t1 - lcdf_diff;
  double t2, t2_sign;

  if (lpdf_l == lpdf_r) t2 = lmath_logaddexps(0, t1, 1, t1_sign);
  else
  {
    t2 = lmath_logaddexpss(t1,
                           2 * (lpdf_diff - lcdf_diff),
                           t1_sign,
                           -1,
                           &t2_sign);

    if ((t2_sign == 1) || (t2 < 0)) t2 = lmath_logaddexps(0,
                                                          t2,
                                                          1,
                                                          t2_sign);
    else t2 = 0;
  }
  *lvar_res = log(var_) + t2;
}

void consolidate_mean_variance(int n, Buffer *buffer, double *mean,
                               double *variance)
{
  int i;

  double laccum = logsumexp(buffer->lw, n);
  double c;

  *mean = 0;

  for (i = 0; i < n; i++)
  {
    c          = exp(buffer->lw[i] - laccum);
    *mean     += c * buffer->m[i];
    *variance += c * buffer->v[i];
  }

  *variance -= (*mean) * (*mean);
}

void mean_variance(void (*f_fl_fll)(double,
                                    void *,
                                    double *,
                                    double *,
                                    double *), void *args,
                   Interval *interval, Normal *normal, Buffer *buffer,
                   double *mean, double *variance)
{
  int i;
  double left, right;

  for (i = 0; i < interval->n - 1; i++)
  {
    left  = i * interval->left;
    right = left + interval->step;
    integrate_window(f_fl_fll, args, left, right, normal,
                     buffer->lw + i, buffer->m + i, buffer->v + i);
  }

  consolidate_mean_variance(interval->n - 1, buffer, mean, variance);
}

//
// void compute_mean_variance(double *N,
//                            double *K,
//                            double *eta,
//                            double *tau,
//                            double *lmom0,
//                            double *mu_res,
//                            double *var_res,
//                            int     N_len,
//                            int     _nintervals,
//                            double *_height,
//                            double *_weight)
// {
//   // global _nintervals, _height, _weight
//   int i;
//   double left, right;
//   double step;
//   double mu, var;
//
//   for (i = 0; i < N_len; i++)
//   {
//     mu  = eta[i] / tau[i];
//     var = 1. / tau[i];
//     meaningful_interval(N[i], K[i], mu, var, &left, &right, _height,
//                         _weight, _nintervals);
//     step = (right - left) / _nintervals;
//     moments(left, step, _nintervals, N[i], K[i], mu, var,
//             &lmom0[i], &mu_res[i], &var_res[i]);
//   }
// }
