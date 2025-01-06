#ifndef LIKNORM_H
#define LIKNORM_H

/** Major liknorm version. */
#define LIKNORM_VERSION_MAJOR 1
/** Minor liknorm version. */
#define LIKNORM_VERSION_MINOR 5
/** Minor liknorm version. */
#define LIKNORM_VERSION_PATCH 8
/** Liknorm version. */
#define LIKNORM_VERSION "1.5.8"

struct LikNormMachine;

/** Create a Machine instance capable of doing numerical integration.
 *
 * @param size Number of integration points. ``500`` points should be enough.
 * ``300`` is usually fine too.
 * @returns Machine instance to perform integration.
 */
struct LikNormMachine *liknorm_create_machine(int size);

/** Perform numerical integration.
 *
 * @param  machine Machine to perform integration.
 * @param log_zeroth Zeroth moment.
 * @param log_mean First moment of the normalized distribution.
 * @param log_variance Variance of the normalized distribution.
 */
void liknorm_integrate(struct LikNormMachine *, double *log_zeroth,
                       double *mean, double *variance);

/** Destroy a Machine instance.
 *
 * @param machine Machine to be destroyed. Always call it before exiting your
 * program, otherwise it will leak memory.
 */
void liknorm_destroy_machine(struct LikNormMachine *);

double liknorm_logprod(struct LikNormMachine *, double x);

/** Bernoulli distribution.
 *
 * It is the discrete probability distribution of a random variable which takes
 * the value ``1`` with probability ``p`` and the value ``0`` with probability
 * ``1 − p``. (Wikipedia.)
 *
 * @param machine Liknorm handler.
 * @param k Number of successes.
 */
void liknorm_set_bernoulli(struct LikNormMachine *, double k);
void liknorm_set_probit(struct LikNormMachine *, double k);

/** Binomial distribution.
 *
 * It is the discrete probability distribution of the number of successes ``k``
 * in a sequence of ``n`` independent experiments. (Wikipedia.) The probability
 * mass function is given by:
 *
 *     Binom(k, n) pᵏ (1 - p)ⁿ⁻ᵏ,
 *
 * for which ``Binom(m, n) = n! / (m! (n - m)!)``.
 *
 * @param machine Liknorm handler.
 * @param k Number of successes.
 * @param n Number of trials.
 */
void liknorm_set_binomial(struct LikNormMachine *, double k, double n);

/** Negative binomial distribution.
 *
 * It is a discrete probability distribution of the number of successes in a
 * sequence of independent and identically distributed Bernoulli trials before a
 * specified (non-random) number of failures (denoted ``r``) occurs.
 * (Wikipedia.) Let ``k`` be the number of successes. The probability mass
 * function is given by:
 *
 *     Binom(k, k + r - 1) pᵏ (1 - p)ʳ,
 *
 * for which ``Binom(m, n) = n! / (m! (n - m)!)``.
 *
 * @param machine Liknorm handler.
 * @param k Number of successes.
 * @param r Number of failures.
 */
void liknorm_set_nbinomial(struct LikNormMachine *, double k, double r);

/** Poisson distribution.
 *
 * It is a discrete probability distribution that expresses the probability of a
 * given number of events occurring in a fixed interval of time or space if
 * these events occur with a known constant rate and independently of the time
 * since the last event. (Wikipedia.) The probability mass function is given by:
 *
 *     λᵏe^{-λ} / k!,
 *
 * for which ``λ`` is the rate of occurrence and ``k`` the number of
 * occurrences.
 *
 * @param machine Liknorm handler.
 * @param k Number of occurrences.
 */
void liknorm_set_poisson(struct LikNormMachine *, double k);

/** Exponential distribution.
 *
 * It is the probability distribution that describes the time between events in
 * a Poisson point process, i.e., a process in which events occur continuously
 * and independently at a constant average rate. (Wikipedia.)
 *
 * @param machine Liknorm handler.
 * @param x Time between events.
 */
void liknorm_set_exponential(struct LikNormMachine *, double x);

/** Gamma distribution.
 *
 * A gamma distribution is a general type of statistical distribution that is
 * related to the beta distribution and arises naturally in processes for which
 * the waiting times between Poisson distributed events are relevant. (Wolfram.)
 *
 * @param machine Liknorm handler.
 * @param x Waiting time.
 * @param a Shape parameter α.
 */
void liknorm_set_gamma(struct LikNormMachine *, double x, double a);

/** Geometric distribution.
 *
 * @param machine Liknorm handler.
 * @param x Number of trials to success.
 */
void liknorm_set_geometric(struct LikNormMachine *, double x);

/** Set the natural parameters of Normal prior.
 *
 * @param machine: Machine to perform integration.
 * @param tau It equals to ``σ⁻²``.
 * @param eta It equals to ``μσ⁻²``.
 */
void liknorm_set_prior(struct LikNormMachine *, double tau, double eta);

#endif
