#include <math.h>
#include "base.h"
#include "moments.h"
#include "ncephes/cprob.h"

void compute_heights(double left, double step, int nint,
                     double N, double K, double mu, double var,
                     double* heights)
{
    int i;
    for (i = 0; i < nint+1; i++)
    {
        heights[i] = log_ulike_prior(left, N, K, mu, var);
        left += step;
    }
}

double compute_weights(double left, double step, int nint,
                       double N, double K, double mu, double var,
                       double* heights,
                       double* w)
{
    double hl, hr, c;
    double r, llc;
    double acum;
    int i;
    for (i = 0; i < nint; i++)
    {
        r = left + step;

        hl = heights[i];
        hr = heights[i+1];
        c = (hr - hl) / step;

        llc = log(fabs(c));
        if (c > 0)
            w[i] = hl - left*c + lmath_logaddexps(r*c, left*c, 1, -1) - llc;
        else
            w[i] = hl - left*c + lmath_logaddexps(r*c, left*c, -1, 1) - llc;

        if (i == 0)
            acum = w[i];
        else
            acum = lmath_logaddexp(acum, w[i]);

        left = r;
    }

    return acum;
}

double shrink_intervals(double* left,
                        double step, int nint,
                        double N, double K, double mu, double var,
                        double* heights, double* w)
{

    compute_heights(left[0], step, nint, N, K, mu, var, heights);
    double r;
    double acum;
    int i;

    acum = compute_weights(left[0], step, nint, N, K, mu, var,
                           heights, w);

    for (i = 0; i < nint; i++)
        w[i] = exp(w[i] - acum);

    i = nint;
    while (w[i-1] <= 1e-256)
        i--;
    r = left[0] + step * i;

    int j = 0;
    while (w[j+1] <= 1e-256 && j < i-1)
        j++;
    left[0] += step * j;

    return r;
}

void moments_array(double* N, double* K,
             double* eta, double* tau,
             double* lmom0, double* mu_res, double* var_res, int N_len,
             int _nintervals, double* _height, double* _weight)
{
    // global _nintervals, _height, _weight
    int i;
    double left, right;
    double step;
    double mu, var;

    for (i = 0; i < N_len; i++)
    {
        mu = eta[i]/tau[i];
        var = 1./tau[i];
        meaningful_interval(N[i], K[i], mu, var, &left, &right, _height,
                            _weight, _nintervals);
        step = (right - left) / _nintervals;
        moments(left, step, _nintervals, N[i], K[i], mu, var,
                &lmom0[i], &mu_res[i], &var_res[i]);
    }
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
