#ifndef POISSON_H
#define POISSON_H

/** Poisson log-partition function.
 *
 * Definition:
 *
 *     b(𝜃) = exp(𝜃)
 */
double poisson_log_partition(const double theta);

/** Log of the first derivative of the Poisson log-partition function.
 *
 * Definition:
 *
 *     log(b'(𝜃)) = 𝜃
 */
double poisson_log_partition_fderivative(const double theta);

/** Log of the derivatives of the Poisson log-partition function.
 *
 * Implements ``b(𝜃)``, ``log(b'(𝜃))``, and:
 *
 *     log(b''(𝜃)) = 𝜃
 *
 */
void poisson_log_partition_derivatives(const double theta, double *b0, double *logb1,
                                       double *logb2);

#endif
