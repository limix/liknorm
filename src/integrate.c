#include "integrate.h"
#include "expfam.h"
#include "logaddexp.h"
#include "machine.h"
#include "normal.h"
#include <assert.h>
#include <float.h>
#include <math.h>

#define LPI2 0.572364942924700081938738094323
#define LNSQRT2 0.346573590279972698624533222755

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

    double beta;
    double alpha;

    double lcdf_a, lcdf_b, lsf_a, lsf_b;
    double lcdf_diff;
    double logpbeta;
    double logpalpha;
    double logp, logp_sign;

    double D;

    double tsxx;
    double tsyy;

    double tsx;
    double tsy;

    double C;
    int Csign;

    double nominator;

    double htau2;

    int sign = htau + tmidiff > 0 ? -1 : +1;

    b += sign * exp(*logA1 + log(fabs(htau + tmidiff)) - log(htau));
    sign = 2 * htau + tmidiff > 0 ? +1 : -1;
    a += sign * tmi *
         exp(*logA1 + log(fabs(2 * htau + tmidiff)) - log(2 * htau));
    assert(isfinite(b));

    assert(isfinite(a));
    assert(isfinite(b));

    beta = (tsii - b) / htau_sqrt;
    alpha = (tsi - b) / htau_sqrt;

    if (alpha + beta >= 0) {
        lsf_a = logcdf(-alpha);
        lsf_b = logcdf(-beta);
        lcdf_diff = lsf_a + log1p(-exp(-lsf_a + lsf_b));
    } else {
        lcdf_a = logcdf(alpha);
        lcdf_b = logcdf(beta);
        lcdf_diff = lcdf_b + log1p(-exp(-lcdf_b + lcdf_a));
    }

    *log_zeroth =
        (a + (b * b) / 2) / htau + LPI2 + LNSQRT2 - log_htau / 2 + lcdf_diff;

    assert(isfinite(*log_zeroth));

    logpbeta = logpdf(beta);
    logpalpha = logpdf(alpha);

    if (logpbeta > logpalpha) {
        logp = logpbeta + log1p(-exp(-logpbeta + logpalpha));
        logp_sign = 1;
    } else {
        logp = logpalpha + log1p(-exp(-logpalpha + logpbeta));
        logp_sign = -1;
    }

    D = logp_sign * exp(logp - lcdf_diff);

    tsxx = log(fabs(b + tsi)) + logp;
    tsyy = log(tstep) + logpbeta;

    tsx = logp_sign * (b + tsi);
    tsy = tstep;

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

    nominator = fmax(b * b + htau - htau_sqrt * C, DBL_EPSILON);

    htau2 = htau * htau;
    assert(htau2 > 0);
    *v = nominator / htau2;

    *u = (htau * (b - htau_sqrt * D)) / htau2;

    assert(isfinite(htau) && htau >= 0);
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
        assert(isfinite(m->diff[i]));
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
        assert(isfinite(m->u[i]));
        assert(isfinite(m->v[i]));
        *mean += m->u[i] * m->diff[i];
        *variance += m->v[i] * m->diff[i];
    }

    *variance = *variance - (*mean) * (*mean);

    assert(isfinite(*variance));
    assert(isfinite(*mean));

    *variance = fmax(*variance, DBL_EPSILON);

    step = (*right - *left) / machine->size;
    *left += ileft * step;
    *right -= (m->size - iright - 1) * step;
}
