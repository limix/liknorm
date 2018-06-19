#ifndef LOGADDEXP_H
#define LOGADDEXP_H

#include "compiler.h"
#include <math.h>

/* Implements log(e^x + e^y).
 */
static inline double logaddexp(const double x, const double y) {
    double m = fmax(x, y);
    return m + log(exp(x - m) + exp(y - m));
}

/* Implements log(e^x - e^y).
 *
 * It assumes that e^x - e^y > 0.
 */
static inline double logsubexp(const double x, const double y) {
    return x + log1p(-exp(y - x));
}

static inline double logaddexp_array(const double *x, const int n,
                                     const double xmax) {
    double total = 0;
    int i;

    for (i = 0; i < n; ++i)
        total += exp(x[i] - xmax);

    return xmax + log(total);
}

#endif
