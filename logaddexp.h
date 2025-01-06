#ifndef LIKNORM_LOGADDEXP_H
#define LIKNORM_LOGADDEXP_H

/* For Windows. */
#define _USE_MATH_DEFINES

#include <float.h>
#include <math.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942 /* log_e 2 */
#endif

/* Computes ㏒ₑ(𝑒ˣ + 𝑒ʸ) in a safe and accurate way.
 *
 * For example, `log(exp(1e3) + exp(-INFINITY))` will likely overflow,
 * while `logaddexp(1e3, -INFINITY)` will return `1e3`.
 */
inline static double logaddexp(double x, double y)
{
    double const tmp = x - y;

    if (x == y) return x + M_LN2;

    if (tmp > 0)
        return x + log1p(exp(-tmp));
    else if (tmp <= 0)
        return y + log1p(exp(tmp));

    return tmp;
}

#endif
