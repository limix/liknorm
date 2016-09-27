#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "constants.h"
#include "normal.h"


typedef struct
{
  double left, right;
  int    n;
  double step;
} Interval;

typedef void log_partition (double  theta,
                            double *b0,
                            double *logb1,
                            double *logb2);

typedef struct
{
        double y;
        double aphi;
        double log_aphi;
        log_partition* lp;
        double left;
        double right;
} ExpFam;

typedef struct
{
        double eta;
        double log_tau;
        double tau;
} Normal;

typedef struct
{
        double *log_zeroth;
        double *u;
        double *v;
        int n;
        double *A0;
        double *logA1;
        double *logA2;
        double *midiff;
        double precision;
} LikNormMachine;


/* Implements log(e^x + e^y).
 */
inline static double logaddexp(double x, double y)
{
    double tmp = x - y;

    if (x == y)
        return x + M_LN2;

    if (tmp > 0)
        return x + log1p(exp(-tmp));
    else if (tmp <= 0)
        return y + log1p(exp(tmp));

    return tmp;
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
void combine_steps(LikNormMachine *machine, double *mean, double *variance);
void shrink_interval(ExpFam *ef, double step, double *left, double *right);
void integrate(LikNormMachine *machine,
               ExpFam         *ef,
               Normal         *normal,
               double         *mean,
               double         *variance);
LikNormMachine* create_liknorm_machine(int n, double precision);
void destroy_liknorm_machine(LikNormMachine *machine);
void binomial_log_partition(double  theta, double *b0, double *logb1,
                            double *logb2);
void bernoulli_log_partition(double  theta, double *b0, double *logb1,
                             double *logb2);

void poisson_log_partition(double  theta, double *b0, double *logb1,
                           double *logb2);

void gamma_log_partition(double theta, double *b0, double *logb1,
                         double *logb2);
void exponential_log_partition(double  theta, double *b0, double *logb1,
                               double *logb2);
void geometric_log_partition(double  theta, double *b0, double *logb1,
                             double *logb2);
log_partition* get_log_partition(const char *name);
void get_interval(const char *name, double *left, double *right);



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

  *v = hmu * hmu + hvar - hstd * copysign(exp(tmp), tmp_sign);
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

void shrink_interval(ExpFam *ef, double step, double *left, double *right)
{
  double b0;
  double limit = 700 * ef->aphi;

  goto left_loop;

  while (*left < *right && fabs(*left * ef->y - b0) > limit)
  {
    *left += step;
left_loop:;
    (*ef->lp)(*left, &b0, 0, 0);

    // printf("fabs(*left * ef->y - b0) %g\n", fabs(*left * ef->y - b0));
  }

  goto right_loop;

  while (*left < *right && fabs(*right * ef->y - b0) > limit)
  {
    *right -= step;
right_loop:;
    (*ef->lp)(*right, &b0, 0, 0);

    // printf("fabs(*right * ef->y - b0) %g\n", fabs(*right * ef->y - b0));
  }
}

void integrate(LikNormMachine *machine,
               ExpFam         *ef,
               Normal         *normal,
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

  // printf("left right %g %g\n", left, right);

  if (left >= ef->right)
  {
    *mean     = 0;
    *variance = 0;
    return;
  }

  double step = (right - left) / machine->n;

  shrink_interval(ef, step, &left, &right);
  step = (right - left) / machine->n;

  // printf("left right %g %g\n", left, right);

  Interval interval;

  interval.left  = left;
  interval.right = right;
  interval.n     = machine->n;
  interval.step  = step;

  for (int i = 0; i < interval.n; ++i)
  {
    double mi = (2 * interval.left + (2 * i + 1) * interval.step) / 2;
    (*ef->lp)(mi, machine->A0 + i, machine->logA1 + i, machine->logA2 + i);
  }

  for (int i = 0; i < interval.n; ++i)
  {
    machine->A0[i]    /= ef->aphi;
    machine->logA1[i] -= ef->log_aphi;
    machine->logA2[i] -= ef->log_aphi;
  }

  for (int i = 0; i < interval.n; ++i)
  {
    double mi = (2 * interval.left + (2 * i + 1) * interval.step) / 2;
    machine->midiff[i] = -mi* exp(machine->logA2[i] - machine->logA1[i]);
  }

  for (int i = 0; i < interval.n; ++i)
  {
    integrate_step(interval.left + interval.step * i,
                   interval.step,
                   ef,
                   normal,
                   machine->log_zeroth + i,
                   machine->u + i,
                   machine->v + i,
                   machine->A0 + i,
                   machine->logA1 + i,
                   machine->logA2 + i,
                   machine->midiff + i);
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
  machine->A0         = malloc(n * sizeof(double));
  machine->logA1      = malloc(n * sizeof(double));
  machine->logA2      = malloc(n * sizeof(double));
  machine->midiff     = malloc(n * sizeof(double));
  machine->precision  = precision;

  return machine;
}

void destroy_liknorm_machine(LikNormMachine *machine)
{
  free(machine->log_zeroth);
  free(machine->u);
  free(machine->v);
  free(machine->A0);
  free(machine->logA1);
  free(machine->logA2);
  free(machine->midiff);
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
  *b0 = theta + log1p(exp(-theta));

  if (logb1 == 0) return;

  *logb1 = theta - *b0;
  *logb2 = theta - 2 * (*b0);
}

void bernoulli_log_partition(double  theta,
                             double *b0,
                             double *logb1,
                             double *logb2)
{
  *b0 = theta + log1p(exp(-theta));

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

  *logb2 = *logb1 = theta;
}

void gamma_log_partition(double  theta,
                         double *b0,
                         double *logb1,
                         double *logb2)
{
  *b0 = -log(-theta);

  if (logb1 == 0) return;

  *logb1 = *b0;
  *logb2 = 2 * (*b0);
}

void exponential_log_partition(double  theta,
                               double *b0,
                               double *logb1,
                               double *logb2)
{
  *b0 = -log(-theta);

  if (logb1 == 0) return;

  *logb1 = *b0;
  *logb2 = 2 * (*b0);
}

void geometric_log_partition(double  theta,
                             double *b0,
                             double *logb1,
                             double *logb2)
{
  *b0 = -log1p(-exp(theta));

  if (logb1 == 0) return;

  *logb1 = theta + *b0;
  *logb2 = theta + 2 * (*b0);
}

log_partition* get_log_partition(const char *name)
{
  if (strcmp(name, "binomial") == 0) return binomial_log_partition;

  if (strcmp(name, "bernoulli") == 0) return bernoulli_log_partition;

  if (strcmp(name, "poisson") == 0) return poisson_log_partition;

  if (strcmp(name, "gamma") == 0) return gamma_log_partition;

  if (strcmp(name, "exponential") == 0) return exponential_log_partition;

  if (strcmp(name, "geometric") == 0) return geometric_log_partition;

  return 0;
}

void get_interval(const char *name, double *left, double *right)
{
  *left  = -DBL_MAX;
  *right = +DBL_MAX;

  if (strcmp(name, "gamma") == 0) *right = -1e-15;

  if (strcmp(name, "exponential") == 0) *right = -1e-15;

  if (strcmp(name, "geometric") == 0) *right = -1e-15;
}
