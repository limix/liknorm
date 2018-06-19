#include "integrate.h"
#include "expfam.h"
#include "logaddexp.h"
#include "machine.h"
#include "normal.h"
#include "sign.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

#define LPI2 0.572364942924700081938738094323
#define LNSQRT2 0.346573590279972698624533222755

void integrate_step(double si, double step, ExpFam *ef, Normal *normal,
                    double *log_zeroth, double *u, double *v, double *A0,
                    double *logA1, double *logA2, double *diff) {

    const double log_htau = logaddexp(normal->log_tau, *logA2);
    const double htau = exp(log_htau);
    const double htau_sqrt = sqrt(htau);

    const double tstep = step * htau;
    const double tsi = si * htau;
    const double tsii = tsi + tstep;

    const double tmi = (tsi + tsii) / 2;
    const double tmidiff = *diff * tmi;

    double a = -(*A0) * htau;
    double b = ef->y / ef->aphi + normal->eta;

    double lcdf_a, lcdf_b, lsf_a, lsf_b;

    double C;
    int Csign;

    double s = -copysign(1, htau + tmidiff);

    b += s * exp(*logA1 + log(fabs(htau + tmidiff)) - log(htau));
    s = copysign(1, 2 * htau + tmidiff);
    a += s * tmi * exp(*logA1 + log(fabs(2 * htau + tmidiff)) - log(2 * htau));

    const double beta = (tsii - b) / htau_sqrt;
    const double alpha = (tsi - b) / htau_sqrt;

    const double lcdf_diff =
        logsubexp(logcdf((alpha < -beta) * (beta + alpha) - alpha),
                  logcdf((alpha < -beta) * (beta + alpha) - beta));

    *log_zeroth =
        (a + (b * b) / 2) / htau + LPI2 + LNSQRT2 - log_htau / 2 + lcdf_diff;

    const double logpbeta = logpdf(beta);
    const double logpalpha = logpdf(alpha);
    const double logp_sign = copysign(1, logpbeta - logpalpha);

    const double logp =
        logsubexp((logpbeta > logpalpha) * (logpbeta - logpalpha) + logpalpha,
                  (logpbeta > logpalpha) * (logpalpha - logpbeta) + logpbeta);

    const double D = logp_sign * exp(logp - lcdf_diff);

    const double tsxx = log(fabs(b + tsi)) + logp;
    const double tsyy = log(tstep) + logpbeta;

    const double tsx = logp_sign * (b + tsi);
    const double tsy = tstep;

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

    const double nominator = fmax(b * b + htau - htau_sqrt * C, DBL_EPSILON);

    const double htau2 = htau * htau;
    assert(htau2 > 0);
    *v = nominator / htau2;

    *u = (htau * (b - htau_sqrt * D)) / htau2;

    assert(*v >= 0);
}

void combine_steps(LikNormMachine *machine, double *log_zeroth, double *mean,
                   double *variance, double *left, double *right) {

    LikNormMachine *m = machine;
    int i, ileft, iright;
    double step;

    double max_log_zeroth = m->log_zeroth[0];
    for (i = 1; i < m->size; ++i)
        max_log_zeroth = fmax(m->log_zeroth[i], max_log_zeroth);

    (*log_zeroth) = logaddexp_array(m->log_zeroth, m->size, max_log_zeroth);

    for (i = 0; i < m->size; ++i) {
        m->diff[i] = exp(m->log_zeroth[i] - *log_zeroth);
    }

    ileft = -1;

    while (m->diff[++ileft] == 0)
        ;

    iright = m->size;

    while (m->diff[--iright] == 0)
        ;
    ++iright;

    assert(ileft < iright);

    *mean = 0;
    *variance = 0;

    for (i = ileft; i < iright; ++i) {
        *mean += m->u[i] * m->diff[i];
        *variance += m->v[i] * m->diff[i];
    }

    *variance = *variance - (*mean) * (*mean);

    *variance = fmax(*variance, DBL_EPSILON);

    step = (*right - *left) / machine->size;
    *left += ileft * step;
    *right -= (m->size - iright - 1) * step;
}
