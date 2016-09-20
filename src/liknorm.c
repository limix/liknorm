#include "liknorm_impl.h"
#include "liknorm.h"
#include "logaddexp.h"
#include "normal.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846


void integrate_step(double     si,
                    double     step,
                    ExpFam     ef,
                    Normal     normal,
                    LogMoments lm)
{
  double sii = si + step;

  double mi  = si / 2 + sii / 2;
  double eta = normal.eta;
  double tau = normal.tau;

  double A0, logA1, logA2, sign;

  (*ef.lp)(mi, ef.lp_data, &A0, &logA1, &logA2, &sign);

  double tmp, tmp_sign;

  tmp = logaddexpss(logA1, logA2, sign * mi, -(mi * mi) / 2, &tmp_sign);
  double a = -A0 + copysign(exp(tmp), tmp_sign);

  tmp = logaddexpss(logA1, logA2, -sign, mi, &tmp_sign);
  double b = ef.Ty + eta + copysign(exp(tmp), tmp_sign);

  double c  = -(tau + exp(logA2)) / 2;
  double bc = b / (2 * c);


  double hmu   = -b / (2 * c);
  double hvar  = 1 / (-2 * c);
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

  // printf(":%f:\n", lcdf_diff);

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

  // printf(" a: %.10f b: %.10f c: %.10f lcdf_diff: %.50f\n", a, b, c,
  // lcdf_diff);
  *(lm.log_zeroth) = a - b * bc / 2 + log(PI) / 2 - log(-c) / 2 + lcdf_diff;

  // printf("  lm.log_zeroth %.10f\n",                        *(lm.log_zeroth));

  double u = hmu - logp_sign * hstd * exp(logp - lcdf_diff);

  *(lm.u) = u;

  // printf("  lm.u      %.10f\n", *(lm.u));

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

  for (i = 0; i < machine->n; i++)
  {
    (*mean)     += machine->u[i] * exp(machine->log_zeroth[i] - total);
    (*variance) += machine->v[i] * exp(machine->log_zeroth[i] - total);
  }
  (*variance) = (*variance) - (*mean) * (*mean);
}

void shrink_interval(ExpFam ef, double step, double *left, double *right)
{
  double A0;

  goto left_loop;

  while (*left < *right && fabs(*left * ef.Ty - A0) > 700)
  {
    *left += step;
left_loop:;
    (*ef.lp)(*left, ef.lp_data, &A0, 0, 0, 0);
  }

  goto right_loop;

  while (*left < *right && fabs(*right * ef.Ty - A0) > 700)
  {
    *right -= step;
right_loop:;
    (*ef.lp)(*right, ef.lp_data, &A0, 0, 0, 0);
  }
}

void integrate(LikNormMachine *machine,
               ExpFam          ef,
               Normal          normal,
               double         *mean,
               double         *variance)
{
  double left  = normal.eta / normal.tau - 10 * sqrt(1 / normal.tau);
  double right = normal.eta / normal.tau + 10 * sqrt(1 / normal.tau);

  double step = (right - left) / machine->n;

  // double before_left  = left;
  // double before_right = right;

  shrink_interval(ef, step, &left, &right);
  step = (right - left) / machine->n;

  // printf("left %.10f -> %.10f\n",    before_left,  left);
  // printf("right %.10f -> %.10f\n\n", before_right, right);

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

  // printf("mean     %.30f\n", *mean);
  // printf("variance %.30f\n", *variance);
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
