#ifndef EXPFAM_H
#define EXPFAM_H

typedef double log_partition(const double theta);

typedef void log_partition_derivatives(const double theta, double *b0, double *logb1,
                                       double *logb2);

typedef double log_partition_fderivative(const double theta);

enum lik_name
{
    liknorm_bernoulli,
    liknorm_binomial,
    liknorm_poisson,
    liknorm_exponential,
    liknorm_gamma,
    liknorm_geometric,
    liknorm_probit
};

/** Exponential family of distributions.
 *
 * We adopt the following representation:
 *
 *     exp{(y𝜃 - b(𝜃))/a(𝜙) + c(y,𝜙)}.
 *
 * Definitions
 * -----------
 *
 * - y is the random variable.
 * - 𝜃 is the canonical parameter
 * - 𝜙 is the nuisance parameter
 * - a(𝜙) TODO
 * - b(𝜃) is the log-partition function
 * - c(y,𝜙)
 *
 * The mean and variance of y are given by:
 *
 *     E[y]   = b'(𝜃)
 *     Var[y] = b''(𝜃)a(𝜙)
 *
 * Of fundamental importance is the natural parameter η. Given a link function g(・),
 * the natural relate to the canonical parameters via the y mean:
 *
 *     η = g(E[y]) = g(b'(𝜃))
 *
 * If g(・) happens to be the canonical link function, we have:
 *
 *     η = 𝜃
 */
struct ExpFam
{
    double y;
    double a;                        /**< a(𝜙) */
    double loga;                     /**< log(a(𝜙)) */
    double c;                        /**< c(y,𝜙) */
    log_partition *lp;               /**< b(𝜃) */
    log_partition_fderivative *lpfd; /**< b'(𝜃) */
    log_partition_derivatives *lpd;  /**< b''(𝜃) */
    double lower_bound;
    double upper_bound;
    enum lik_name name;
};

#endif
