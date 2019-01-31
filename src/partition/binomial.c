#include "partition/binomial.h"
#include "compiler.h"

#include <float.h>
#include <math.h>

/** Binomial log-partition function.
 *
 * Definition:
 *
 *     b(𝜃) = log(1 + exp(𝜃)).
 */
double binomial_log_partition(const double theta)
{
    return theta < log(DBL_MAX) ? log1p(exp(theta)) : theta;
}

/** First derivative of the Binomial log-partition function.
 *
 * Definition
 * ----------
 *
 * - b'(𝜃) = exp(𝜃) / (1 + exp(𝜃))
 * - log(b'(𝜃)) = 𝜃 - log(1 + exp(𝜃))
 */
double binomial_log_partition_fderivative(const double theta)
{
    return theta < log(DBL_MAX) ? theta - log1p(exp(theta)) : 0;
}

/** Zeroth, first, and second derivatives of the Binomial log-partition function.
 *
 * Definition
 * ----------
 *
 * - b(𝜃) = log(1 + exp(𝜃))
 * - log(b'(𝜃)) = 𝜃 - log(1+exp(𝜃))
 * - log(b''(𝜃)) = 𝜃 - 2log(1+exp(𝜃))
 */
void binomial_log_partition_derivatives(const double theta, double *b0, double *logb1,
                                        double *logb2)
{
    if (theta < log(DBL_MAX)) {
        double log1p_ = log1p(exp(theta));
        *b0 = log1p_;
        *logb1 = theta - log1p_;
        *logb2 = theta - 2 * log1p_;
    } else {
        *b0 = theta;
        *logb1 = 0;
        *logb2 = -theta;
    }
}
