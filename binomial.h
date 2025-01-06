#ifndef BINOMIAL_H
#define BINOMIAL_H

/** Binomial log-partition function.
 *
 * Definition:
 *
 *     b(𝜃) = log(1 + exp(𝜃)).
 */
double liknorm_binomial_log_partition(const double theta);

/** First derivative of the Binomial log-partition function.
 *
 * Definition:
 *
 *     log(b'(𝜃)) = 𝜃 - log(1 + exp(𝜃))
 */
double liknorm_binomial_log_partition_fderivative(const double theta);

/** Zeroth, first, and second derivatives of the Binomial log-partition
 * function.
 *
 * Implements ``b(𝜃)``, ``log(b'(𝜃))``, and:
 *
 *     log(b''(𝜃)) = 𝜃 - 2log(1 + exp(𝜃))
 */
void liknorm_binomial_log_partition_derivatives(const double theta, double *b0,
                                                double *logb1, double *logb2);

#endif
