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

void meaningful_interval(double N, double K, double mu, double var,
                         double* left, double* right, double* heights,
                         double* w, int w_len)
{
    get_extrema(N, K, mu, var, left, right);

    int nint = w_len;
    double step = (right[0] - left[0]) / nint;

    right[0] = shrink_intervals(left, step, nint, N, K, mu, var, heights, w);

    step = (right[0] - left[0]) / nint;
    right[0] = shrink_intervals(left, step, nint, N, K, mu, var, heights, w);
}

// # import numpy as np
// # cimport numpy as np
// from cpython cimport array
// import array
//
// cdef int _nintervals
// cdef double[:] _height
// cdef double[:] _weight
//
// cpdef init(int nintervals):
//     global _nintervals, _height, _weight
//     _nintervals = nintervals
//
//     _height = array.clone(array.array('d', []), nintervals+1, zero=False)
//     _weight = array.clone(array.array('d', []), nintervals, zero=False)

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
