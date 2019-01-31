#define LIKNORM_API_EXPORTS

#include "liknorm.h"
#include "compiler.h"
#include "integrate.h"
#include "interval.h"
#include "machine.h"
#include "partition/partition.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

LIKNORM_API struct LikNormMachine *liknorm_create_machine(int size)
{
    struct LikNormMachine *machine = malloc(sizeof(struct LikNormMachine));

    machine->size = size;
    machine->log_zeroth = malloc(size * sizeof(double));
    machine->u = malloc(size * sizeof(double));
    machine->v = malloc(size * sizeof(double));
    machine->A0 = malloc(size * sizeof(double));
    machine->logA1 = malloc(size * sizeof(double));
    machine->logA2 = malloc(size * sizeof(double));
    machine->diff = malloc(size * sizeof(double));

    return machine;
}

void _liknorm_integrate(struct LikNormMachine *machine, double *log_zeroth,
                        double *mean, double *variance, double *left, double *right)
{
    struct ExpFam *ef = &(machine->ef);
    struct Normal *normal = &(machine->normal);
    int i;

    double step = (*right - *left) / machine->size;
    double *A0 = machine->A0;
    double *logA1 = machine->logA1;
    double *logA2 = machine->logA2;
    double *diff = machine->diff;
    double *u;
    double *v;
    double *mlog_zeroth;
    const double pi = 3.14159265358979323846;

    for (i = 0; i < machine->size; ++i)
        (*ef->lpd)(*left + step * i + step / 2, A0 + i, logA1 + i, logA2 + i);

    for (i = 0; i < machine->size; ++i) {
        A0[i] /= ef->a;
        logA1[i] -= ef->loga;
        logA2[i] -= ef->loga;
        diff[i] = -exp(logA2[i] - logA1[i]);
    }

    u = machine->u;
    v = machine->v;
    mlog_zeroth = machine->log_zeroth;

    for (i = 0; i < machine->size; ++i) {
        integrate_step(*left + step * i, step, ef, normal, mlog_zeroth++, u++, v++,
                       A0++, logA1++, logA2++, diff++);
    }

    combine_steps(machine, log_zeroth, mean, variance, left, right);

    *log_zeroth += machine->ef.c;
    *log_zeroth -= log((2 * pi) / normal->tau) / 2;
    *log_zeroth -= (normal->eta * normal->eta) / (2 * normal->tau);
}

void liknorm_integrate_probit(double y, double tau, double eta, double *log_zeroth,
                              double *mean, double *variance)
{
    double c, b, denom, logpdfc, logcdfc, logdiff;
    double tau1 = tau + 1;
    double d = sqrt(tau) / sqrt(tau1);
    y = 2 * y - 1;
    c = (sqrt(tau) * y * eta / sqrt(tau1)) / tau;
    logcdfc = logcdf(c);
    logpdfc = logpdf(c);
    *log_zeroth = logcdfc;
    logdiff = exp(logpdfc - logcdfc);
    b = logdiff + c;
    denom = 1 - b * logdiff / (tau + 1);
    *variance = denom / tau;
    *mean = (eta + y * logdiff * d) / (1 - b * logdiff / tau1);
    *mean = *mean * (*variance);
}

LIKNORM_API void liknorm_integrate(struct LikNormMachine *machine, double *log_zeroth,
                                   double *mean, double *variance)
{

    double left, right;
    struct ExpFam *ef = &(machine->ef);
    struct Normal *normal = &(machine->normal);
    double ileft;
    double iright;

    if (ef->name == liknorm_probit) {
        liknorm_integrate_probit(ef->y, normal->tau, normal->eta, log_zeroth, mean,
                                 variance);
        return;
    }

    find_interval(ef, normal, &left, &right);

    do {
        ileft = left;
        iright = right;
        _liknorm_integrate(machine, log_zeroth, mean, variance, &left, &right);
    } while ((right - left) / (iright - ileft) < 0.9);
}

LIKNORM_API void liknorm_destroy_machine(struct LikNormMachine *machine)
{
    free(machine->log_zeroth);
    free(machine->u);
    free(machine->v);
    free(machine->A0);
    free(machine->logA1);
    free(machine->logA2);
    free(machine->diff);
    free(machine);
}

LIKNORM_API double liknorm_logprod(struct LikNormMachine *machine, double x)
{

    const double PI = 3.141592653589793238462643383279502884;
    double a = machine->ef.a;
    double b = machine->ef.lp(x);
    double c = machine->ef.c;
    double y = machine->ef.y;
    double left = (y * x - b) / a + c;
    double tau = machine->normal.tau;
    double eta = machine->normal.eta;
    double right = -(log(2 * PI) - machine->normal.log_tau) / 2;
    right += -x * x * tau / 2 + eta * x - eta * eta / (2 * tau);
    return left + right;
}

/** Binomial distribution.
 *
 * We are assuming the use of the canonical link function in this description.
 *
 * Definitons
 * ----------
 *
 * - k is the success (1) of failure (0).
 * - y       = k
 * - 𝜙       = 1
 * - a(𝜙)    = 𝜙 = 1
 * - b(𝜃)    = log(1 + exp(𝜃))
 * - c(y,𝜙) = log(binom(𝜙, 𝜙y)) = log(binom(n, k))
 *
 * The support is therefore y ϵ {0/n, 1/n, ..., n/n}.
 */
LIKNORM_API void liknorm_set_bernoulli(struct LikNormMachine *machine, double k)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_bernoulli;
    m->ef.y = k;
    m->ef.a = 1;
    m->ef.loga = 0;
    m->ef.c = 0;
    m->ef.lp = bernoulli_log_partition;
    m->ef.lpfd = bernoulli_log_partition_fderivative;
    m->ef.lpd = bernoulli_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = +DBL_MAX;
}

LIKNORM_API void liknorm_set_probit(struct LikNormMachine *machine, double k)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_probit;
    m->ef.y = k;
}

double logbinom(double k, double n)
{
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

/** Binomial distribution.
 *
 * We are assuming the use of the canonical link function (logit).
 *
 * Link functions
 * --------------
 *
 * Let p = E[y]. We have
 *
 * - canonical(p)     = logit(p)    = log(p/(1-p))
 * - canonical_inv(η) = logistic(η) = 1 / (1 + exp(-η))
 *
 * Definitons
 * ----------
 *
 * - n is the number of trials
 * - k is the number of successes
 * - y       = k/n
 * - 𝜙       = n
 * - a(𝜙)   = 1/𝜙 = 1/n
 * - b(𝜃)    = log(1 + exp(𝜃))
 * - b'(𝜃)   = exp(𝜃) / (1 + exp(𝜃))
 * - b''(𝜃)  = exp(𝜃) / (1+exp(𝜃))²
 * - c(y,𝜙) = log(binom(𝜙, 𝜙y)) = log(binom(n, k))
 *
 * The support is therefore y ϵ {0/n, 1/n, ..., n/n}.
 */
LIKNORM_API void liknorm_set_binomial(struct LikNormMachine *machine, double k,
                                      double n)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_binomial;
    m->ef.y = k / n;
    m->ef.a = 1 / n;
    m->ef.loga = -log(n);
    m->ef.c = logbinom(k, n);
    m->ef.lp = binomial_log_partition;
    m->ef.lpfd = binomial_log_partition_fderivative;
    m->ef.lpd = binomial_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = +DBL_MAX;
}

/** Negative binomial distribution.
 *
 * We are assume the canonical link function.
 *
 * @param k Number of successes.
 * @param r Number of failures.
 *
 *
 * Notes
 * -----
 *
 * The exponential family functions are:
 *
 *     a(𝜙) = 1/𝜙
 *     b(θ) = -log(1 - exp(θ))
 *     b'(θ) = exp(θ) / (1 - exp(θ))
 *     b''(θ) = exp(θ) / (1 - exp(θ))²
 *     c(y,𝜙) = log(binom(y𝜙 + 𝜙 - 1, y𝜙))
 *
 *
 * y = k/r
 * 𝜙 = r
 * If p is the probability of TODO, we have
 *
 *     θ = log(1 - p)
 *
 * The support is therefore y ϵ {0/n, 1/n, ..., n/n}.
 * Link functions
 * --------------
 *
 * Let 𝜇 = E[y]. We have
 *
 * - canonical(𝜇)     = log(𝜇/(r+𝜇))
 * - canonical_inv(η) =
 */
LIKNORM_API void liknorm_set_nbinomial(struct LikNormMachine *machine, double k,
                                       double r)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_binomial;
    m->ef.y = k;
    m->ef.a = 1;
    m->ef.loga = 0;
    m->ef.c = logbinom(k, k + r - 1);
    m->ef.lp = nbinomial_log_partition;
    m->ef.lpfd = nbinomial_log_partition_fderivative;
    m->ef.lpd = nbinomial_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = +DBL_MAX;
}

static inline double logfactorial(double k) { return lgamma(k + 1); }

/** Poisson distribution.
 *
 * We are assuming the use of the canonical link function.
 *
 * Link functions
 * --------------
 *
 * Let λ = E[y]. We have
 *
 * - canonical(λ)     = log(λ)
 * - canonical_inv(η) = exp(η)
 *
 * Definitons
 * ----------
 *
 * - k is the number of occurrences
 * - y        = k
 * - 𝜙       = 1
 * - a(𝜙)    = 𝜙 = 1
 * - b(𝜃)     = exp(𝜃)
 * - b'(𝜃)    = exp(𝜃)
 * - b'(𝜃)    = exp(𝜃)
 * - c(y,𝜙)  = -log(y!) = -log(k!)
 *
 * The support is therefore y ϵ {0, 1, 2, ...}.
 */
LIKNORM_API void liknorm_set_poisson(struct LikNormMachine *machine, double k)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_poisson;
    m->ef.y = k;
    m->ef.a = 1;
    m->ef.loga = 0;
    m->ef.c = -logfactorial(k);
    m->ef.lp = poisson_log_partition;
    m->ef.lpfd = poisson_log_partition_fderivative;
    m->ef.lpd = poisson_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = +DBL_MAX;
}

LIKNORM_API void liknorm_set_exponential(struct LikNormMachine *machine, double x)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_exponential;
    m->ef.y = x;
    m->ef.a = 1;
    m->ef.loga = 0;
    m->ef.c = 0;
    m->ef.lp = exponential_log_partition;
    m->ef.lpfd = exponential_log_partition_fderivative;
    m->ef.lpd = exponential_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = -DBL_EPSILON;
}

LIKNORM_API void liknorm_set_gamma(struct LikNormMachine *machine, double x, double a)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_gamma;
    m->ef.y = x;
    m->ef.a = 1 / a;
    m->ef.loga = -log(a);
    m->ef.c = 0;
    m->ef.lp = gamma_log_partition;
    m->ef.lpfd = gamma_log_partition_fderivative;
    m->ef.lpd = gamma_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = -DBL_EPSILON;
}

LIKNORM_API void liknorm_set_geometric(struct LikNormMachine *machine, double x)
{
    struct LikNormMachine *m = machine;
    m->ef.name = liknorm_geometric;
    m->ef.y = x;
    m->ef.a = 1;
    m->ef.loga = 0;
    m->ef.c = 0;
    m->ef.lp = geometric_log_partition;
    m->ef.lpfd = geometric_log_partition_fderivative;
    m->ef.lpd = geometric_log_partition_derivatives;
    m->ef.lower_bound = -DBL_MAX;
    m->ef.upper_bound = -DBL_EPSILON;
}

LIKNORM_API void liknorm_set_prior(struct LikNormMachine *machine, double tau,
                                   double eta)
{
    const double tau_min = 2 * sqrt(DBL_EPSILON);
    tau = fmax(tau, tau_min);
    machine->normal.eta = eta;
    machine->normal.tau = tau;
    machine->normal.log_tau = log(tau);
}
