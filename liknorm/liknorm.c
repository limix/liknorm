#include "liknorm.h"
#include "machine.h"
#include <stdlib.h>

enum lik_name {
  bernoulli,
  binomial,
  poisson,
  exponential,
  gamma,
  geometric
};

static const char lik_name_str[] = {
  "bernoulli",
  "binomial",
  "poisson",
  "exponential",
  "gamma",
  "geometric"
};

LikNormMachine *liknorm_create_machine(int size) {
  LikNormMachine *machine = malloc(sizeof(LikNormMachine));

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

void liknorm_destroy_machine(LikNormMachine *machine) {
  free(machine->log_zeroth);
  free(machine->u);
  free(machine->v);
  free(machine->A0);
  free(machine->logA1);
  free(machine->logA2);
  free(machine->diff);
  free(machine);
}

void liknorm_set_bernoulli(LikNormMachine *machine, double k) {
  // set_expfam(machine, BERNOULLI, k, 1, 0, "bernoulli");
  LikNormMachine *m = machine;
  m->ef.name = bernoulli;
  m->ef.y = k;
  m->ef.aphi = 1;
  m->ef.log_aphi = 0;
  m->ef.c = 0;
  m->ef.lp = get_log_partition(name);
  m->ef.lp0 = get_log_partition0(name);
  m->ef.lp1 = get_log_partition1(name);
  get_likelihood_interval(name, &m->ef.lower_bound, &m->ef.upper_bound);
}

void liknorm_set_binomial(LikNormMachine *machine, double k, double n) {
  set_expfam(machine, BINOMIAL, k / n, 1 / n, logbinom(k, n),
             "binomial");
}

void liknorm_set_poisson(LikNormMachine *machine, double k) {
  set_expfam(machine, POISSON, k, 1, -logfactorial(k),
             "poisson");
}
void liknorm_set_exponential(LikNormMachine *machine, double x) {
  set_expfam(machine, EXPONENTIAL, x, 1, 0, "exponential");
}
void liknorm_set_gamma(LikNormMachine *machine, double x, double a) {
  set_expfam(machine, GAMMA, x, 1 / a, 0, "gamma");
}
void liknorm_set_geometric(LikNormMachine *machine, double x) {
  set_expfam(machine, GEOMETRIC, x, 1, 0, "geometric");
}


void liknorm_set_prior(LikNormMachine *machine, double tau, double eta) {
  assert(tau > 0);
  tau = fmax(tau, LIK_SQRT_EPSILON * 2);
  machine->normal.eta = eta;
  machine->normal.tau = tau;
  machine->normal.log_tau = log(tau);
}
