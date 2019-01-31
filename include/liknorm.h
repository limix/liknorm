#ifndef LIKNORM_H
#define LIKNORM_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef LIKNORM_API_EXPORTS
#define LIKNORM_API __declspec(dllexport)
#else
#define LIKNORM_API __declspec(dllimport)
#endif
#else
#define LIKNORM_API
#endif

struct LikNormMachine;

LIKNORM_API struct LikNormMachine *liknorm_create_machine(int size);
LIKNORM_API void liknorm_integrate(struct LikNormMachine *machine,
                                   double *log_zeroth, double *mean,
                                   double *variance);
LIKNORM_API void liknorm_destroy_machine(struct LikNormMachine *machine);

LIKNORM_API double liknorm_logprod(struct LikNormMachine *machine, double x);

LIKNORM_API void liknorm_set_bernoulli(struct LikNormMachine *machine,
                                       double k);
LIKNORM_API void liknorm_set_probit(struct LikNormMachine *machine, double k);
LIKNORM_API void liknorm_set_binomial(struct LikNormMachine *machine, double k,
                                      double n);
LIKNORM_API void liknorm_set_poisson(struct LikNormMachine *machine, double k);
LIKNORM_API void liknorm_set_exponential(struct LikNormMachine *machine,
                                         double x);
LIKNORM_API void liknorm_set_gamma(struct LikNormMachine *machine, double x,
                                   double a);
LIKNORM_API void liknorm_set_geometric(struct LikNormMachine *machine,
                                       double x);
LIKNORM_API void liknorm_set_prior(struct LikNormMachine *machine, double tau,
                                   double eta);

#ifdef __cplusplus
}
#endif

#endif
