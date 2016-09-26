#include "liknorm_impl.h"
#include "liknorm.h"
#include "logaddexp.h"
#include "normal.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define PI 3.14159265358979323846

// log(PI) / 2
#define LPI2 0.572364942924700081938738094323


void integrate_step(double     si,
                    double     step,
                    ExpFam     ef,
                    Normal     normal,
                    LogMoments lm)
{
  double sii = si + step;

  double mi = (si + sii) / 2;

  double b0, logb1, logb2;

  (*ef.lp)(mi, &b0, &logb1, &logb2);

  double A0    = b0 / ef.aphi;
  double logA1 = logb1 - ef.log_aphi;
  double logA2 = logb2 - ef.log_aphi;

  // int mi_pos = mi > 0;

  double tmp, tmp_sign;

  // tmp = logaddexpss(logA1, logA2, mi, -(mi * mi) / 2, &tmp_sign);
  // double a = -A0;

  // A1 + (-mi/2) A2
  bool   positive = mi < 0 || log(mi / 2) + logA2 < logA1;
  double diff     = logA1 - logA2;
  double a        = -A0;

  if (positive) a += mi * exp(logA1 + log1p((-mi / 2) * exp(-diff)));
  else a -= mi * exp(logA2 + log(mi / 2) + log1p((-2 / mi) * exp(+diff)));

  // a = -A0 + copysign(exp(tmp), tmp_sign);

  tmp = logaddexpss(logA1, logA2, -1, mi, &tmp_sign);
  double Ty = ef.y / ef.aphi;
  double b  = Ty + normal.eta + copysign(exp(tmp), tmp_sign);

  double c = -(normal.tau + exp(logA2)) / 2;


  double hvar  = -1 / (2 * c);
  double hmu   = b * hvar;
  double hstd  = sqrt(hvar);
  double beta  = (sii - hmu) / hstd;
  double alpha = (si - hmu) / hstd;

  double lcdf_a, lcdf_b, lsf_a, lsf_b;
  double lcdf_diff;

  if (alpha + beta >= 0)
  {
    lsf_a = logcdf(-alpha);
    lsf_b = logcdf(-beta);

    lcdf_diff = logaddexps(lsf_a, lsf_b, 1.0, -1.0);
  } else {
    lcdf_a = logcdf(alpha);
    lcdf_b = logcdf(beta);

    lcdf_diff = logaddexps(lcdf_b, lcdf_a, 1.0, -1.0);
  }

  double logpbeta  = logpdf(beta);
  double logpalpha = logpdf(alpha);

  double logp, logp_sign;

  if (fabs(beta) < fabs(alpha))
  {
    logp      = logaddexps(logpbeta, logpalpha, +1, -1);
    logp_sign = +1;
  }
  else {
    logp      = logaddexps(logpbeta, logpalpha, -1, +1);
    logp_sign = -1;
  }

  *(lm.log_zeroth) = a + (b * hmu) / 2 + LPI2 - log(-c) / 2 + lcdf_diff;

  *(lm.u) = hmu - logp_sign * hstd * exp(logp - lcdf_diff);

  tmp  = logaddexpss(logpbeta, logpalpha, hmu + sii, -hmu - si, &tmp_sign);
  tmp -= lcdf_diff;

  *(lm.v) = hmu * hmu + hvar - hstd * copysign(exp(tmp), tmp_sign);
}

void combine_steps(LikNormMachine *machine, double *mean, double *variance)
{
  (*mean)     = 0;
  (*variance) = 0;

  double total = logaddexp(machine->log_zeroth[0], machine->log_zeroth[1]);
  int    i;

  for (i = 2; i < machine->n; i++) total = logaddexp(total,
                                                     machine->log_zeroth[i]);

  double diff;

  for (i = 0; i < machine->n; i++)
  {
    diff         = exp(machine->log_zeroth[i] - total);
    (*mean)     += machine->u[i] * diff;
    (*variance) += machine->v[i] * diff;
  }
  (*variance) = (*variance) - (*mean) * (*mean);
}

void shrink_interval(ExpFam ef, double step, double *left, double *right)
{
  double b0;

  goto left_loop;

  while (*left < *right && fabs(*left * ef.y / ef.aphi - b0 / ef.aphi) > 700)
  {
    *left += step;
left_loop:;
    (*ef.lp)(*left, &b0, 0, 0);
  }

  goto right_loop;

  while (*left < *right && fabs(*right * ef.y / ef.aphi - b0 / ef.aphi) > 700)
  {
    *right -= step;
right_loop:;
    (*ef.lp)(*right, &b0, 0, 0);
  }
}

void integrate(LikNormMachine *machine,
               ExpFam          ef,
               Normal          normal,
               double         *mean,
               double         *variance)
{
  double left = normal.eta / normal.tau - 10 * sqrt(1 / normal.tau);

  left = fmax(left, ef.left);

  double right = normal.eta / normal.tau + 10 * sqrt(1 / normal.tau);
  right = fmin(right, ef.right);

  if (left >= ef.right)
  {
    *mean     = 0;
    *variance = 0;
    return;
  }

  double step = (right - left) / machine->n;

  shrink_interval(ef, step, &left, &right);
  step = (right - left) / machine->n;

  Interval   interval;
  LogMoments lm;

  interval.left  = left;
  interval.right = right;
  interval.n     = machine->n;
  interval.step  = step;

  double si;

  for (int i = 0; i < interval.n; i++)
  {
    si            = interval.left + interval.step * i;
    lm.log_zeroth = machine->log_zeroth + i;
    lm.u          = machine->u + i;
    lm.v          = machine->v + i;

    integrate_step(si, interval.step, ef, normal, lm);
  }

  combine_steps(machine, mean, variance);
}

LikNormMachine* create_liknorm_machine(int n, double precision)
{
  LikNormMachine *machine = malloc(sizeof(LikNormMachine));

  machine->n          = n;
  machine->log_zeroth = malloc(n * sizeof(double));
  machine->u          = malloc(n * sizeof(double));
  machine->v          = malloc(n * sizeof(double));
  machine->precision  = precision;

  return machine;
}

void destroy_liknorm_machine(LikNormMachine *machine)
{
  free(machine->log_zeroth);
  free(machine->u);
  free(machine->v);
  free(machine);
}

// \phi = N
// a(\phi) = 1/\phi
// b(\theta) = log(1 + e^\theta)
// b'(\theta) = e^\theta / (1 + e^\theta)
// b''(\theta) = e^\theta / (1 + e^\theta)^2
void binomial_log_partition(double  theta,
                            double *b0,
                            double *logb1,
                            double *logb2)
{
  *b0 = logaddexp(0, theta);

  if (logb1 == 0) return;

  *logb1 = theta - *b0;
  *logb2 = theta - 2 * (*b0);
}

void bernoulli_log_partition(double  theta,
                             double *b0,
                             double *logb1,
                             double *logb2)
{
  *b0 = logaddexp(0, theta);

  if (logb1 == 0) return;

  *logb1 = theta - *b0;
  *logb2 = theta - 2 * (*b0);
}

void poisson_log_partition(double  theta,
                           double *b0,
                           double *logb1,
                           double *logb2)
{
  *b0 = exp(theta);

  if (logb1 == 0) return;

  if (logb1 != 0) *logb1 = theta;

  if (logb2 != 0) *logb2 = theta;
}

void gamma_log_partition(double  theta,
                         double *b0,
                         double *logb1,
                         double *logb2)
{
  *b0 = log(-1 / theta);

  if (logb1 == 0) return;

  *logb1 = -log(-theta);
  *logb2 = -2 * log(fabs(theta));
}

void exponential_log_partition(double  theta,
                               double *b0,
                               double *logb1,
                               double *logb2)
{
  *b0 = log(-1 / theta);

  if (logb1 == 0) return;

  *logb1 = -log(-theta);
  *logb2 = -2 * log(fabs(theta));
}

void geometric_log_partition(double  theta,
                             double *b0,
                             double *logb1,
                             double *logb2)
{
  *b0 = -logaddexps(0, theta, 1, -1);

  if (logb1 == 0) return;

  *logb1 = theta + *b0;
  *logb2 = theta + 2 * (*b0);
}

log_partition* get_log_partition(char *name)
{
  if (strcmp(name, "binomial") == 0) return binomial_log_partition;

  if (strcmp(name, "bernoulli") == 0) return bernoulli_log_partition;

  if (strcmp(name, "poisson") == 0) return poisson_log_partition;

  if (strcmp(name, "gamma") == 0) return gamma_log_partition;

  if (strcmp(name, "exponential") == 0) return exponential_log_partition;

  if (strcmp(name, "geometric") == 0) return geometric_log_partition;

  return 0;
}

void get_interval(char *name, double *left, double *right)
{
  *left  = -DBL_MAX;
  *right = +DBL_MAX;

  if (strcmp(name, "gamma") == 0) *right = -1e-15;

  if (strcmp(name, "exponential") == 0) *right = -1e-15;

  if (strcmp(name, "geometric") == 0) *right = -1e-15;
}
