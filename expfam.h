#ifndef LIKNORM_EXPFAM_H
#define LIKNORM_EXPFAM_H

typedef double liknorm_log_partition(const double theta);

typedef void liknorm_log_partition_derivatives(const double theta, double *b0,
                                       double *logb1, double *logb2);

typedef double liknorm_log_partition_fderivative(const double theta);

enum liknorm_likelihood
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
 * We adopt the following representation::
 *
 *     f(y; θ, 𝜙) = exp{(yθ - b(θ))/a(𝜙) + c(y,𝜙)},
 *
 * for which::
 *
 *     y     : random variable value.
 *     θ     : canonical parameter.
 *     𝜙     : nuisance parameter.
 *     a(𝜙)  :
 *     b(θ)  : log-partition function.
 *     c(y,𝜙): normaliser.
 *
 * The mean and variance are given by::
 *
 *     E[y]   = b'(θ)
 *     Var[y] = b''(θ)a(𝜙)
 *
 * In order to define a generalised linear mixed model (GLMM) we use the
 * so-called natural parameter ``η``. Given a link function ``g(.)``, the
 * natural parameter relates to the canonical parameter as follows::
 *
 *     η = g(E[y]) = g(b'(θ)).
 *
 * Every member of the exponential family has a canonical link function, which
 * greatly simplifies the relationship::
 *
 *     η = θ
 */
struct ExpFam
{
    double y;                        /**< Random variable value */
    double a;                        /**< ``a(𝜙)`` */
    double loga;                     /**< ``log(a(𝜙))`` */
    double c;                        /**< ``c(y,𝜙)`` */
    liknorm_log_partition *lp;               /**< ``b(θ)`` */
    liknorm_log_partition_fderivative *lpfd; /**< ``log(b'(θ))`` */
    liknorm_log_partition_derivatives *lpd;  /**< ``log(b''(θ))`` */
    double lower_bound;
    double upper_bound;
    enum liknorm_likelihood name;
};

#endif
