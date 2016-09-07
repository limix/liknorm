#include <stdio.h>
#include <math.h>
#include "base.h"

void get_extrema(double N, double K, double mu, double var,
                 double* left, double* right)
{
  *left = mu - 9*sqrt(var);
  *right = mu + 9*sqrt(var);

  double c = get_mode(N, K);
  double ileft = c - 10.0;
  double iright = c + 10.0;
  if ((*left < ileft) && (iright < *right))
  {
    *left = ileft;
    *right = iright;
  }
  else
  {
    *left = fmin(ileft, *left);
    *right = fmax(iright, *right);
  }
}

void update_moms(double lmom0_int, double lmu_int, double lmu_int_sign,
                 double lvar_int, double* lmom0_out, double* lmu_out,
                 double* lmu_out_sign, double* lvar_out)
{
  double sign;
  *lmom0_out = lmath_logaddexp(*lmom0_out, lmom0_int);

  *lmu_out = lmath_logaddexpss(*lmu_out, lmu_int + lmom0_int, *lmu_out_sign,
                         lmu_int_sign, &sign);
  *lmu_out_sign = sign;

  *lvar_out = lmath_logaddexp(*lvar_out, lvar_int + lmom0_int);
  *lvar_out = lmath_logaddexp(*lvar_out, 2*lmu_int + lmom0_int);
}

void tail_integral(double x, int left, double mu, double var,
                   double* lmom0, double* lmu_res, double* lmu_res_sign,
                   double* lvar_res)
{
  double sig = sqrt(var);
  double lsig = log(sig);
  double m = (x - mu)/sig;
  double lpdf = lmath_normal_logpdf(m);

  if (left)
    *lmom0 = lmath_normal_logsf(m);
  else
    *lmom0 = lmath_normal_logcdf(m);

  double sign;
  if (left)
    sign = +1.0;
  else
    sign = -1.0;

  *lmu_res = mu + sign * exp(lpdf - *lmom0 + lsig);
  if (*lmu_res >= 0.)
    *lmu_res_sign = +1.0;
  else
  {
    *lmu_res_sign = -1.0;
    *lmu_res *= -1.0;
  }
  *lmu_res = log(*lmu_res);

  *lvar_res = var + sign * sig * (mu + x) * exp(lpdf - *lmom0)
              - exp(2**lmu_res) + mu*mu;
  *lvar_res = log(*lvar_res);
}

void f_fl_fll(double x, double N, double K,
              double* f_, double* fl_, double* fll_)
{
  double lpdf = lmath_normal_logpdf(x);
  double lcdf = lmath_normal_logcdf(x);
  double lsf = lmath_normal_logsf(x);
  double a = 0, b = 0;
  double lK = 0, lNK = 0;
  double sign, lr;

  *f_ = K * lcdf + (N-K) * lsf;

  double l1 = lpdf - lcdf;
  double l2 = lpdf - lsf;

  if (K > 0)
  {
      lK = log(K);
      a = lK + l1;
  }

  if (K < N)
  {
    lNK = log(N-K);
    b = lNK + l2;
  }

  if (K > 0 && K < N)
  {
    if (a == b)
      *fl_ = 0.0;
    else
    {
      lr = lmath_logaddexpss(a, b, 1.0, -1.0, &sign);
      *fl_ = sign * exp(lr);
    }
  }
  else
  {
    if (K == 0)
      *fl_ = -exp(b);
    else
      *fl_ = exp(a);
  }

  if (K == 0)
  {
    a = exp(lNK + l2);
    b = exp(lNK + 2*l2);
  }
  else if (K == N)
  {
    a = - exp(lK + l1);
    b = + exp(lK + 2*l1);
  }
  else
  {
    a = - exp(lK + l1)   + exp(lNK + l2);
    b = + exp(lK + 2*l1) + exp(lNK + 2*l2);
  }
  *fll_ = x * a - b;
}

inline double normal_logpdf2(double x, double mu, double var)
{
    return lmath_normal_logpdf((x - mu) / sqrt(var)) - log(var)/2.0;
}

double log_ulike_prior(double x, double N, double K, double mu, double var)
{
  double la = lmath_normal_logcdf(x);
  double lb = lmath_normal_logsf(x);
  return (K * la + (N-K) * lb + normal_logpdf2(x, mu, var));
}

void integrate_window(double l_, double r_,
                      double N, double K,
                      double mu, double var,
                      double* lmom0, double* lmu_res, double* lmu_res_sign,
                      double* lvar_res)
{
  double x0 = l_ + (r_-l_)/2.0;
  double a, b, c;
  // printf("int 1\n"); fflush(stdout);
  f_fl_fll(x0, N, K, &a, &b, &c);

  double alpha = a - b*x0 + c*x0*x0/2.0 - mu*mu/(2*var);
  double beta = b - c*x0 + mu/var;
  double gamma = c/2.0 - 1.0/(2*var);


  double mu_ = -beta/(2*gamma);
  double var_ = -1./(2*gamma);
  double sig_ = sqrt(var_);


  double lcoef = -beta*beta/(4*gamma) + alpha;
  lcoef -= log(-2 * gamma * var) / 2.0;

  double ml = (l_ - mu_)/sig_;
  double mr = (r_ - mu_)/sig_;
  // printf("int 2 mr %.10f ml %.10f\n", mr, ml); fflush(stdout);

  double lpdf_r = lmath_normal_logpdf(mr);
  double lpdf_l = lmath_normal_logpdf(ml);
  double lcdf_l, lcdf_r;
  double lcdf_diff, lpdf_diff;
  double lsf_l, lsf_r;
  // printf("int 2.1 lpdf_r %.10f lpdf_l %.10f\n", lpdf_r, lpdf_l); fflush(stdout);

  if (ml + mr >= 0.0)
  {
    lsf_l = lmath_normal_logsf(ml);
    lsf_r = lmath_normal_logsf(mr);
    // printf("lsf_l lsf_r %.10f %.10f\n", lsf_l, lsf_r); fflush(stdout);
    lcdf_diff = lmath_logaddexps(lsf_l, lsf_r, 1.0, -1.0);
    // printf("lcdf_diff %.10f\n", lcdf_diff); fflush(stdout);
  }
  else
  {
    lcdf_l = lmath_normal_logcdf(ml);
    lcdf_r = lmath_normal_logcdf(mr);
    lcdf_diff = lmath_logaddexps(lcdf_r, lcdf_l, 1.0, -1.0);
  }
  // printf("int 3\n"); fflush(stdout);

  *lmom0 = lcoef + lcdf_diff;
  double lpdf_diff_sign;

  if (lpdf_l == lpdf_r)
  {
    *lmu_res = log(fabs(mu_));
    if (mu_ >= 0) *lmu_res_sign = +1.0;
    else *lmu_res_sign = -1.0;
  }
  else
  {
    lpdf_diff = lmath_logaddexpss(lpdf_l, lpdf_r, +1.0, -1.0, &lpdf_diff_sign);
    *lmu_res = mu_ + sig_ * exp(lpdf_diff - lcdf_diff) * lpdf_diff_sign;
    if (*lmu_res >= 0) *lmu_res_sign = +1.0;
    else *lmu_res_sign = -1.0;
    *lmu_res = log(fabs(*lmu_res));
  }
  // printf("int 4\n"); fflush(stdout);
  double t1_sign;
  double t1 = lmath_logaddexpss(lpdf_l, lpdf_r, ml, -mr, &t1_sign);
  t1 = t1 - lcdf_diff;
  double t2, t2_sign;

  // printf("int 4.1 t1 %.10f\n", t1); fflush(stdout);

  if (lpdf_l == lpdf_r)
    t2 = lmath_logaddexps(0.0, t1, 1.0, t1_sign);
  else
  {
    // printf("int 4.11 lpdf_diff lcdf_diff %.10f %.10f %.10f\n", lpdf_diff, lcdf_diff, 2*(lpdf_diff - lcdf_diff)); fflush(stdout);
    // printf("int 4.12 %.40f\n", t1 - 2*(lpdf_diff - lcdf_diff)); fflush(stdout);
    // printf("int 4.13 %.3f\n", t1_sign); fflush(stdout);
    t2 = lmath_logaddexpss(t1, 2*(lpdf_diff - lcdf_diff), t1_sign, -1.0, &t2_sign);
    // printf("int 4.2 t2 %.10f %.1f\n", t2, t2_sign); fflush(stdout);
    if (t2_sign == 1.0 || t2 < 0.0)
      t2 = lmath_logaddexps(0.0, t2, 1.0, t2_sign);
    else
      t2 = 0.0;
    // printf("int 4.3 t2 %.10f\n", t2); fflush(stdout);
  }
  // printf("int 5\n"); fflush(stdout);
  *lvar_res = log(var_) + t2;
}

void moments(double left, double step, int nints,
             double N, double K, double mu, double var,
             double* lmom0_res, double* mu_res, double* var_res)
{
  int i;
  double l_ = left;
  double r_;
  double lmom0_int;
  double lmu_int, lmu_int_sign;
  double lvar_int;
  double lmom0_final;
  double lmu_final, lmu_final_sign;
  double lvar_final;

  r_ = l_ + step;
  // printf("ponto 1\n"); fflush(stdout);
  integrate_window(l_, r_, N, K, mu, var,
         &lmom0_final, &lmu_final, &lmu_final_sign, &lvar_final);
  // printf("%.10f ", lmom0_final)
  // printf("ponto 2\n"); fflush(stdout);
  lvar_final = lmath_logaddexp(lvar_final, 2*lmu_final);
  lvar_final += lmom0_final;
  lmu_final += lmom0_final;
  l_ = r_;

  // printf("ponto 3\n"); fflush(stdout);
  for (i = 0; i < nints-1; i++)
  {
    r_ = l_ + step;
    // printf("ponto 3.1 l_ %.10f r_ %.10f\n", l_, r_); fflush(stdout);
    integrate_window(l_, r_, N, K, mu, var, &lmom0_int, &lmu_int,
                     &lmu_int_sign, &lvar_int);
    // printf("ponto 3.2 %.30f\n", lvar_int); fflush(stdout);
    update_moms(lmom0_int, lmu_int, lmu_int_sign, lvar_int,
                &lmom0_final, &lmu_final, &lmu_final_sign, &lvar_final);
    l_ = r_;
  }
  // printf("%.10f ", lmom0_final)
  // printf("ponto 4\n"); fflush(stdout);
  if (N == K)
  {
    // printf("N==K\n");
    tail_integral(left + nints*step, 1, mu, var,
              &lmom0_int, &lmu_int, &lmu_int_sign, &lvar_int);
    update_moms(lmom0_int, lmu_int, lmu_int_sign, lvar_int,
            &lmom0_final, &lmu_final, &lmu_final_sign, &lvar_final);
  }
  else if (K == 0)
  {
    // printf("K==0\n");
    tail_integral(left, 0, mu, var,
              &lmom0_int, &lmu_int, &lmu_int_sign, &lvar_int);
    update_moms(lmom0_int, lmu_int, lmu_int_sign, lvar_int,
            &lmom0_final, &lmu_final, &lmu_final_sign, &lvar_final);
  }
  lmu_final -= lmom0_final;
  lvar_final -= lmom0_final;
  lvar_final = lmath_logaddexps(lvar_final, 2*lmu_final, 1., -1.);
  lmom0_final += lmath_logbinom(N, K);
  // printf("lmath_logbinom(N, K): %.10f ", lmath_logbinom(N, K))
  // printf("%.10f ", lmom0_final)

  *lmom0_res = lmom0_final;
  *mu_res = lmu_final_sign*exp(lmu_final);
  *var_res = exp(lvar_final);
}
