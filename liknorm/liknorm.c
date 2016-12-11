#ifndef LIKNORM_H
#define LIKNORM_H

// #include <stdio.h>
// #include <math.h>
// #include <stdlib.h>
// #include <string.h>
// #include <stdbool.h>
// #include <assert.h>

// #include "constants.h"
// #include "normal.h"
// #include "definitions.h"
// #include "compiler.h"
// #include "optimizer.h"

void liknorm_integrate(LikNormMachine *machine,
                       ExpFam         *ef,
                       Normal         *normal,
                       double         *log_zeroth,
                       double         *mean,
                       double         *variance);
LikNormMachine* liknorm_create_machine(int n);
void            liknorm_destroy_machine(LikNormMachine *machine);

/* -log(DBL_TRUE_MIN) */
#define DBL_LOG_MAX +744.4400719213812180896638892590999603271484375

/* sqrt(DBL_EPSILON) */
#define SQRT_DBL_EPSILON 1.490116119384765625e-08

// #define GOLDEN 2 / (1 + sqrt(5.0))

/* Implements log(e^x + e^y).
 */
inline static double logaddexp(double x, double y)
{
  double tmp = x - y;

  if (LIKNORM_UNLIKELY(x == y)) return x + M_LN2;

  if (tmp > 0) return x + log1p(exp(-tmp));
  else if (tmp <= 0) return y + log1p(exp(tmp));

  return tmp;
}

double logaddexp_array(double *x, int n, double xmax)
{
  double total = 0;

  for (int i = 0; i < n; ++i) total += exp(x[i] - xmax);

  return xmax + log(total);
}

void integrate_step(double  si,
                    double  step,
                    ExpFam *ef,
                    Normal *normal,
                    double *log_zeroth,
                    double *u,
                    double *v,
                    double *A0,
                    double *logA1,
                    double *logA2,
                    double *midiff);
void combine_steps(LikNormMachine *machine,
                   double          max_log_zeroth,
                   double         *log_zeroth,
                   double         *mean,
                   double         *variance);
void smart_shrink_interval(int     n,
                           ExpFam *ef,
                           Normal *normal,
                           double *left,
                           double *right);

void integrate_step(double  si,
                    double  step,
                    ExpFam *ef,
                    Normal *normal,
                    double *log_zeroth,
                    double *u,
                    double *v,
                    double *A0,
                    double *logA1,
                    double *logA2,
                    double *midiff)
{
  double sii = si + step;

  double mi = (si + sii) / 2;

  double tmp, tmp_sign;

  double a     = -(*A0);
  double Ty    = ef->y / ef->aphi;
  double b     = Ty + normal->eta;
  double logmi = log(fabs(mi));

  double falta = *logA1 - logmi - *logA2 + M_LN2;

  if ((mi < 0) || (falta > M_LN2))
  {
    a += mi * exp(*logA1 + log1p(*midiff / 2));
    b -= exp(*logA1 + log1p(*midiff));
  } else
  {
    if (falta > 0)
    {
      a += mi * exp(*logA1 + log1p(*midiff / 2));
      b += exp(*logA2 + logmi + log1p(1 / *midiff));
    } else {
      a -= mi * exp(*logA2 + logmi - M_LN2 + log1p(2 / *midiff));
      b += exp(*logA2 + logmi + log1p(1 / *midiff));
    }
  }

  double hvar = exp(-logaddexp(normal->log_tau, *logA2));

  double hmu   = b * hvar;
  double hstd  = sqrt(hvar);
  double beta  = (sii - hmu) / hstd;
  double alpha = (si - hmu) / hstd;

  double lcdf_a, lcdf_b, lsf_a, lsf_b;
  double lcdf_diff;

  if (alpha + beta >= 0)
  {
    lsf_a     = logcdf(-alpha);
    lsf_b     = logcdf(-beta);
    lcdf_diff = lsf_a + log1p(-exp(-lsf_a + lsf_b));
  } else {
    lcdf_a    = logcdf(alpha);
    lcdf_b    = logcdf(beta);
    lcdf_diff = lcdf_b + log1p(-exp(-lcdf_b + lcdf_a));
  }

  double logpbeta  = logpdf(beta);
  double logpalpha = logpdf(alpha);

  double logp, logp_sign;

  if (logpbeta > logpalpha)
  {
    logp      = logpbeta + log1p(-exp(-logpbeta + logpalpha));
    logp_sign = 1;
  } else {
    logp      = logpalpha + log1p(-exp(-logpalpha + logpbeta));
    logp_sign = -1;
  }

  /* log(pi)/2 */
  double logpi_2 = 0.572364942924700081938738094323;
  *log_zeroth = a + (b * hmu) / 2 + logpi_2 + log(M_SQRT2 * hstd) + lcdf_diff;

  *u = hmu - logp_sign * hstd * exp(logp - lcdf_diff);

  double k = hmu + si;

  double sxx = log(fabs(logp_sign * k)) + logp;
  double syy = log(step) + logpbeta;
  double sx  = logp_sign * k;
  double sy  = step;

  if (sxx > syy)
  {
    if (sx >= 0)
    {
      tmp      = syy + log1p((sx / sy) * exp(logp - logpbeta));
      tmp_sign = +1;
    } else {
      tmp      = sxx + log1p((sy / sx) * exp(-logp + logpbeta));
      tmp_sign = -1;
    }
  } else {
    tmp      = syy + log1p((sx / sy) * exp(logp - logpbeta));
    tmp_sign = +1;
  }

  tmp -= lcdf_diff;

  *v = hvar + (hmu * hmu - hstd * copysign(exp(tmp), tmp_sign));

  assert(isfinite(hvar));
  assert(*v >= 0);
}

void combine_steps(LikNormMachine *machine,
                   double          max_log_zeroth,
                   double         *log_zeroth,
                   double         *mean,
                   double         *variance)
{
  double n          = machine->n;
  LikNormMachine *m = machine;

  (*log_zeroth) = logaddexp_array(m->log_zeroth, m->n, max_log_zeroth);

  for (int i = 0; i < n; ++i)
  {
    m->diff[i] = exp(m->log_zeroth[i] - *log_zeroth);
    assert(isfinite(m->diff[i]));
  }

  int left = -1;

  while (m->diff[++left] == 0) ;

  int right = n;

  while (m->diff[--right] == 0) ;
  ++right;

  assert(left < right);

  *mean     = 0;
  *variance = 0;

  for (int i = left; i < right; ++i)
  {
    assert(isfinite(m->u[i]));
    assert(isfinite(m->v[i]));
    *mean     += m->u[i] * m->diff[i];
    *variance += m->v[i] * m->diff[i];
  }

  if ((right - left) / ((double)n) < 0.10) *variance = 1e-8;
  else *variance = *variance - (*mean) * (*mean);

  assert(isfinite(*variance));
  assert(isfinite(*mean));
}

void fprintf_expfam(FILE *stream, const ExpFam *ef)
{
  fprintf(stream, "ExpFam:\n");
  fprintf(stream, "  name    : %s\n", ef->name);
  fprintf(stream, "  y       : %g\n", ef->y);
  fprintf(stream, "  aphi    : %g\n", ef->aphi);
  fprintf(stream, "  log_aphi: %g\n", ef->log_aphi);
}

void fprintf_normal(FILE *stream, const Normal *normal)
{
  fprintf(stream, "Normal:\n");
  fprintf(stream, "  tau: %g\n", normal->tau);
  fprintf(stream, "  eta: %g\n", normal->eta);
}

static inline double g_function(ExpFam *ef, Normal *normal, double x)
{
  return x * (ef->y / ef->aphi + normal->eta) - (normal->eta * x * x) / 2 -
         ef->lp0(x) / ef->aphi;
}

static inline double g_func_base(double x, void *args)
{
  void  **args_  = args;
  ExpFam *ef     = args_[0];
  Normal *normal = args_[1];

  return g_function(ef, normal, x);
}

static inline double g_derivative(ExpFam *ef, Normal *normal, double x)
{
  return ef->y / ef->aphi + normal->eta - x * normal->tau
         - exp(ef->lp1(x)) / ef->aphi;
}

void smart_shrink_interval(int n, ExpFam *ef,
                           Normal *normal,
                           double *left,
                           double *right)
{
  double gleft  = g_function(ef, normal, *left);
  double gright = g_function(ef, normal, *right);
  double dleft  = g_derivative(ef, normal, *left);
  double dright = g_derivative(ef, normal, *right);
  double gmax;
  double step = (*right - *left) / n;

  void *args[] = { ef, normal };

  if (dleft * dright > 0)
  {
    gmax = dleft > 0 ? gright : gleft;
  } else {
    gmax =
      g_function(ef, normal,
                 zero(*left, *right, 1e-5 * fmax(fabs(gleft), fabs(gright)),
                      &g_func_base, args));
  }

  while (gmax -
         g_function(ef, normal,
                    *left + step) > DBL_LOG_MAX &&
         *left < *right) *left += step;

  *left = fmin(*left, *right);

  while (gmax -
         g_function(ef, normal,
                    *right - step) > DBL_LOG_MAX &&
         *left < *right) *right -= step;
  *right = fmax(*right, *left);
}

void liknorm_integrate(LikNormMachine *machine,
                       ExpFam         *ef,
                       Normal         *normal,
                       double         *log_zeroth,
                       double         *mean,
                       double         *variance)
{
  const double times = 7;
  double std         = sqrt(1 / normal->tau);
  double mu          = normal->eta / normal->tau;
  double left        = mu - times * std;
  double right       = mu + times * std;

  if (left >= ef->right)
  {
    *mean     = ef->right;
    *variance = 0;
    return;
  }

  if (right <= ef->left)
  {
    *mean     = ef->left;
    *variance = 0;
    return;
  }

  left  = fmax(left, ef->left);
  right = fmin(right, ef->right);


  // printf("BEFORE: left right: %g %g\n", left, right);
  smart_shrink_interval(machine->n, ef, normal, &left, &right);

  // printf("AFTER : left right: %.50f %.50f\n", left, right);
  // printf("right-left: %g\n\n",                right - left);

  if (right - left <= SQRT_DBL_EPSILON)
  {
    *log_zeroth = -DBL_MAX;
    *mean       = left;
    *variance   = 1e-8;
    return;
  }

  double step = (right - left) / machine->n;

  for (int i = 0; i < machine->n; ++i)
  {
    double mi = (2 * left + (2 * i + 1) * step) / 2;
    (*ef->lp)(mi, machine->A0 + i, machine->logA1 + i, machine->logA2 + i);
    machine->A0[i]    /= ef->aphi;
    machine->logA1[i] -= ef->log_aphi;
    machine->logA2[i] -= ef->log_aphi;
    machine->diff[i]   = -mi* exp(machine->logA2[i] - machine->logA1[i]);
  }

  double max_log_zeroth = -DBL_MAX;

  for (int i = 0; i < machine->n; ++i)
  {
    integrate_step(left + step * i, step, ef, normal, machine->log_zeroth + i,
                   machine->u + i,
                   machine->v + i,
                   machine->A0 + i,
                   machine->logA1 + i,
                   machine->logA2 + i,
                   machine->diff + i);
    max_log_zeroth = fmax(max_log_zeroth, machine->log_zeroth[i]);
  }

  combine_steps(machine, max_log_zeroth, log_zeroth, mean, variance);
  *variance = fmax(*variance, 1e-8);
}

LikNormMachine* liknorm_create_machine(int n)
{
  LikNormMachine *machine = malloc(sizeof(LikNormMachine));

  machine->n          = n;
  machine->log_zeroth = malloc(n * sizeof(double));
  machine->u          = malloc(n * sizeof(double));
  machine->v          = malloc(n * sizeof(double));
  machine->A0         = malloc(n * sizeof(double));
  machine->logA1      = malloc(n * sizeof(double));
  machine->logA2      = malloc(n * sizeof(double));
  machine->diff       = malloc(n * sizeof(double));

  return machine;
}

void liknorm_destroy_machine(LikNormMachine *machine)
{
  free(machine->log_zeroth);
  free(machine->u);
  free(machine->v);
  free(machine->A0);
  free(machine->logA1);
  free(machine->logA2);
  free(machine->diff);
  free(machine);
}

#include "log_partitions.h"

#endif /* end of include guard: LIKNORM_H */
