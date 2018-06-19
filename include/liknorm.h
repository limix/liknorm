#ifndef LIKNORM_H
#define LIKNORM_H

struct LikNormMachine;

struct LikNormMachine *liknorm_create_machine(int size);
void liknorm_integrate(struct LikNormMachine *machine, double *log_zeroth,
                       double *mean, double *variance);
void liknorm_destroy_machine(struct LikNormMachine *machine);

void liknorm_set_bernoulli(struct LikNormMachine *machine, double k);
void liknorm_set_probit(struct LikNormMachine *machine, double k);
void liknorm_set_binomial(struct LikNormMachine *machine, double k, double n);
void liknorm_set_poisson(struct LikNormMachine *machine, double k);
void liknorm_set_exponential(struct LikNormMachine *machine, double x);
void liknorm_set_gamma(struct LikNormMachine *machine, double x, double a);
void liknorm_set_geometric(struct LikNormMachine *machine, double x);

void liknorm_set_prior(struct LikNormMachine *machine, double tau, double eta);

#endif
