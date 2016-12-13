#include "expfam.h"
#include "integrate.h"
#include "logaddexp.h"
#include "machine.h"
#include "normal.h"
#include <assert.h>
#include <float.h>
#include <math.h>

void integrate_step(double si, double step, ExpFam *ef, Normal *normal,
                    double *log_zeroth, double *u, double *v, double *A0,
                    double *logA1, double *logA2, double *diff) {

  double log_htau = logaddexp(normal->log_tau, *logA2);
  double htau = exp(log_htau);
  double htau_sqrt = sqrt(htau);

  double tstep = step * htau;
  double tsi = si * htau;
  double tsii = tsi + tstep;

  double tmi = (tsi + tsii) / 2;
  double tmidiff = *diff * tmi;

  double a = -(*A0) * htau;
  double b = ef->y / ef->aphi + normal->eta;

  int sign = htau + tmidiff > 0 ? -1 : +1;
  b += sign * exp(*logA1 + log(fabs(htau + tmidiff)) - log(htau));
  sign = 2 * htau + tmidiff > 0 ? +1 : -1;
  a += sign * tmi * exp(*logA1 + log(fabs(2 * htau + tmidiff)) - log(2 * htau));
  assert(isfinite(b));

  assert(isfinite(a));
  assert(isfinite(b));
  const double heta = b;

  const double beta = (tsii - heta) / htau_sqrt;
  const double alpha = (tsi - heta) / htau_sqrt;

  double lcdf_a, lcdf_b, lsf_a, lsf_b;
  double lcdf_diff;

  if (alpha + beta >= 0) {
    lsf_a = logcdf(-alpha);
    lsf_b = logcdf(-beta);
    lcdf_diff = lsf_a + log1p(-exp(-lsf_a + lsf_b));
  } else {
    lcdf_a = logcdf(alpha);
    lcdf_b = logcdf(beta);
    lcdf_diff = lcdf_b + log1p(-exp(-lcdf_b + lcdf_a));
  }

  /* log(pi)/2 */
  static const double lpi2 = 0.572364942924700081938738094323;
  static const double lnsqrt2 = 0.346573590279972698624533222755;
  *log_zeroth = (a + (heta * heta) / 2) / htau + lpi2 + lnsqrt2 -
                log_htau / 2 + lcdf_diff;

  assert(isfinite(*log_zeroth));

  const double logpbeta = logpdf(beta);
  const double logpalpha = logpdf(alpha);

  double logp, logp_sign;

  if (logpbeta > logpalpha) {
    logp = logpbeta + log1p(-exp(-logpbeta + logpalpha));
    logp_sign = 1;
  } else {
    logp = logpalpha + log1p(-exp(-logpalpha + logpbeta));
    logp_sign = -1;
  }

  const double D = logp_sign * exp(logp - lcdf_diff);

  const double tsxx = log(fabs(heta + tsi)) + logp;
  const double tsyy = log(tstep) + logpbeta;

  const double tsx = logp_sign * (heta + tsi);
  const double tsy = tstep;

  double C;
  int Csign;

  if (tsxx > tsyy) {
    if (tsx >= 0) {
      C = tsyy + log1p((tsx / tsy) * exp(logp - logpbeta));
      Csign = +1;
    } else {
      C = tsxx + log1p((tsy / tsx) * exp(-logp + logpbeta));
      Csign = -1;
    }
  } else {
    C = tsyy + log1p((tsx / tsy) * exp(logp - logpbeta));
    Csign = +1;
  }

  C = Csign * exp(C - lcdf_diff);

  const double nominator =
      fmax(heta * heta + htau - htau_sqrt * C, DBL_EPSILON);

  const double htau2 = htau * htau;
  assert(htau2 > 0);
  *v = nominator / htau2;

  *u = (htau * (heta - htau_sqrt * D)) / htau2;

  assert(isfinite(htau) && htau >= 0);
  assert(*v >= 0);
}

void combine_steps(LikNormMachine *machine, double *log_zeroth, double *mean,
                   double *variance) {

  LikNormMachine *m = machine;

  double max_log_zeroth = m->log_zeroth[0];
  for (int i = 1; i < m->size; ++i)
    max_log_zeroth = fmax(m->log_zeroth[i], max_log_zeroth);

  (*log_zeroth) = logaddexp_array(m->log_zeroth, m->size, max_log_zeroth);

  for (int i = 0; i < m->size; ++i) {
    m->diff[i] = exp(m->log_zeroth[i] - *log_zeroth);
    assert(isfinite(m->diff[i]));
  }

  int left = -1;

  while (m->diff[++left] == 0)
    ;

  int right = m->size;

  while (m->diff[--right] == 0)
    ;
  ++right;

  assert(left < right);

  *mean = 0;
  *variance = 0;

  for (int i = left; i < right; ++i) {
    assert(isfinite(m->u[i]));
    assert(isfinite(m->v[i]));
    *mean += m->u[i] * m->diff[i];
    *variance += m->v[i] * m->diff[i];
  }

  *variance = *variance - (*mean) * (*mean);

  assert(isfinite(*variance));
  assert(isfinite(*mean));
}
