#ifndef LIKNORM_H
#define LIKNORM_H

LikNormMachine* liknorm_create_machine(int n);
void liknorm_integrate(LikNormMachine *machine,
                       double         *log_zeroth,
                       double         *mean,
                       double         *variance);
void liknorm_destroy_machine(LikNormMachine *machine);

void liknorm_set_bernoulli(LikNormMachine *machine, double k);

/* Binomial likelihood.
   k: number of successes
   n: number of trials
*/
void liknorm_set_binomial(LikNormMachine *machine, double k, double n);
void liknorm_set_poisson(LikNormMachine *machine, double k);
void liknorm_set_exponential(LikNormMachine *machine, double x);
void liknorm_set_gamma(LikNormMachine *machine, double x, double a);
void liknorm_set_geometric(LikNormMachine *machine, double x);

void liknorm_set_prior(LikNormMachine *machine, double tau, double eta);

#endif
