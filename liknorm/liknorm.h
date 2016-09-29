#ifndef LIKNORM_H
#define LIKNORM_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "constants.h"
#include "normal.h"
#include "definitions.h"
#include "compiler.h"

/* ========================== Interface ========================== */
void liknorm_integrate(LikNormMachine *machine,
                       ExpFam         *ef,
                       Normal         *normal,
                       double         *log_zeroth,
                       double         *mean,
                       double         *variance);
LikNormMachine* liknorm_create_machine(int n);
void            liknorm_destroy_machine(LikNormMachine *machine);


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

inline static double logaddexps(double x, double y, double sx, double sy)
{
  double tmp = x - y;

  double sxx = log(fabs(sx)) + x;
  double syy = log(fabs(sy)) + y;

  if (sxx == syy)
  {
    if (sx * sy > 0) return sxx + M_LN2;

    return -DBL_MAX;
  }

  if ((sx > 0) && (sy > 0))
  {
    if (tmp > 0) return sxx + log1p((sy / sx) * exp(-tmp));
    else if (tmp <= 0) return syy + log1p((sx / sy) * exp(tmp));
  }
  else if (sx > 0) return sxx + log1p((sy / sx) * exp(-tmp));
  else return syy + log1p((sx / sy) * exp(tmp));

  return tmp;
}

static inline double max_array(double *arr, size_t n)
{
  double v = -DBL_MAX;

  for (size_t i = 0; i < n; ++i) {
    v = fmax(arr[i], v);
  }
  return v;
}

double logaddexps_array(double *x, double *sx, int n)
{
  double total = 0;

  double xmax = max_array(x, n);

  for (int i = 0; i < n; ++i)
  {
    assert(xmax >= x[i]);
    total += sx[i] * exp(x[i] - xmax);
  }

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
                   double          max_log_v,
                   double         *log_zeroth,
                   double         *mean,
                   double         *variance);
int shrink_interval(ExpFam *ef,
                    double  step,
                    double *left,
                    double *right);

/* Returns log(|c|) and c/|c|, for c = sx * e^x + sy * e^y.
 */
inline static double logaddexpss(double x, double y, double sx, double sy,
                                 double *sign)
{
  double sxx = log(fabs(sx)) + x;
  double syy = log(fabs(sy)) + y;

  if (sxx == syy)
  {
    if (sx * sy > 0)
    {
      if (sx > 0) *sign = +1.0;
      else *sign = -1.0;
      return sxx + M_LN2;
    }
    else
    {
      *sign = 1.0;
      return -DBL_MAX;
    }
  }

  if (sxx > syy)
  {
    if (sx >= 0.0) *sign = +1.0;
    else *sign = -1.0;
  }
  else
  {
    if (sy >= 0.0) *sign = +1.0;
    else *sign = -1.0;
  }

  sx *= *sign;
  sy *= *sign;
  return logaddexps(x, y, sx, sy);
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

  *log_zeroth = a + (b * hmu) / 2 + LIK_LPI2 + log(M_SQRT2 * hstd) + lcdf_diff;

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


  // if (ahmu > 1)
  // {
  //   double tmp_c;
  //   tmp_c = hstd * copysign(exp(tmp - log(ahmu)),
  //                           tmp_sign);
  //   *v = hvar + ahmu * (ahmu - tmp_c);
  // } else {
  //   *v = hvar + (hmu * hmu - hstd * copysign(exp(tmp), tmp_sign));
  // }

  double ahmu = fabs(hmu);

  // if ((fabs(tmp) > 1) && (ahmu > 0))
  // {
  //   double porra;
  //   double r = logaddexpss(2 * log(ahmu), tmp, 1, -hstd * tmp_sign, &porra);
  //   *v = logaddexps(log(hvar), r, 1, porra);
  // } else {
  //   *v = log(hvar + (hmu * hmu - hstd * copysign(exp(o), tmp_sign)));
  //
  //   // *v = hvar + (hmu * hmu - hstd * copysign(exp(tmp), tmp_sign));
  // }
  *v = log(hvar + (hmu * hmu - hstd * copysign(exp(tmp), tmp_sign)));


  // *v = log(hvar + (hmu * hmu - hstd * copysign(exp(tmp), tmp_sign)));
  // assert(isfinite(*v));
  assert(isfinite(hvar));
  assert(hvar >= 0);

  // assert((hmu * hmu - hstd * copysign(exp(tmp), tmp_sign)) <= 0);
  assert(hvar + (hmu * hmu - hstd * copysign(exp(tmp), tmp_sign)) > 0);

  // *v = hmu * hmu + hvar - hstd * copysign(exp(tmp), tmp_sign);
}

void combine_steps(LikNormMachine *machine,
                   double          max_log_zeroth,
                   double          max_log_v,
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

  max_log_v = max_array(m->v + left, right - left);

  assert(left < right);

  if ((left > 0) || (right < n))
  {
    printf("Debug: left right perc: %d %d %g%%\n",
           left,
           right,
           (100. * (right - left)) / n);
  }

  (*mean) = 0;

  for (int i = left; i < right; ++i)
  {
    assert(isfinite(m->u[i]));
    (*mean) += m->u[i] * m->diff[i];
  }

  (*variance) = logaddexps_array(m->v + left,
                                 m->diff + left,
                                 right - left);

  if (!isfinite(*variance))
  {
    printf("Non finite variance!!!!!\n");
  }

  double amean = fabs(*mean);

  if (amean > 1)
  {
    (*variance) = exp(logaddexps((*variance), 2 * log(fabs(*mean)), 1, -1));
  }
  else {
    (*variance) = exp(*variance) - (*mean) * (*mean);
  }
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

int shrink_interval(ExpFam *ef, double step, double *left, double *right)
{
  double b0;
  double limit = 7000;

  goto left_loop;

  while (*left < *right && fabs(*left * ef->y - b0) > limit)
  {
    *left += step;
left_loop:;
    b0 = (*ef->lp0)(*left);
  }

  goto right_loop;

  while (*left < *right && fabs(*right * ef->y - b0) > limit)
  {
    *right -= step;
right_loop:;
    b0 = (*ef->lp0)(*right);
  }

  if (*left >= *right)
  {
    fprintf(stderr, "Invalid shrinked interval: [%g, %g].\n", *left, *right);
    return 1;
  }
  return 0;
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

  left = fmax(left, ef->left);

  double right = mu + times * std;
  right = fmin(right, ef->right);


  if (left >= ef->right)
  {
    *mean     = 0;
    *variance = 0;
    return;
  }

  double step = (right - left) / machine->n;

  int err = shrink_interval(ef, step, &left, &right);

  if (err)
  {
    fprintf_normal(stderr, normal);
    fprintf_expfam(stderr, ef);
    exit(EXIT_FAILURE);
  }

  step = (right - left) / machine->n;

  for (int i = 0; i < machine->n; ++i)
  {
    double mi = (2 * left + (2 * i + 1) * step) / 2;
    (*ef->lp)(mi, machine->A0 + i, machine->logA1 + i, machine->logA2 + i);
  }

  for (int i = 0; i < machine->n; ++i)
  {
    machine->A0[i]    /= ef->aphi;
    machine->logA1[i] -= ef->log_aphi;
    machine->logA2[i] -= ef->log_aphi;
  }

  for (int i = 0; i < machine->n; ++i)
  {
    double mi = (2 * left + (2 * i + 1) * step) / 2;
    machine->diff[i] = -mi* exp(machine->logA2[i] - machine->logA1[i]);
  }

  double max_log_zeroth = -DBL_MAX;
  double max_log_v      = -DBL_MAX;

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
    max_log_v      = fmax(max_log_v, machine->v[i]);
  }

  combine_steps(machine, max_log_zeroth, max_log_v, log_zeroth, mean, variance);
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
