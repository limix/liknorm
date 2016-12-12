#include "liknorm.h"
#include "machine.h"
#include "partition/partition.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

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
  LikNormMachine *m = machine;
  m->ef.name = liknorm_bernoulli;
  m->ef.y = k;
  m->ef.aphi = 1;
  m->ef.log_aphi = 0;
  m->ef.c = 0;
  m->ef.lpd = bernoulli_log_partition_derivatives;
  m->ef.lpfd = bernoulli_log_partition_fderivative;
  m->ef.lower_bound = -DBL_MAX;
  m->ef.upper_bound = +DBL_MAX;
}

double logbinom(double k, double n) {
  return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

void liknorm_set_binomial(LikNormMachine *machine, double k, double n) {
  LikNormMachine *m = machine;
  m->ef.name = liknorm_binomial;
  m->ef.y = k / n;
  m->ef.aphi = 1 / n;
  m->ef.log_aphi = -log(n);
  m->ef.c = logbinom(k, n);
  m->ef.lpd = binomial_log_partition_derivatives;
  m->ef.lpfd = binomial_log_partition_fderivative;
  m->ef.lower_bound = -DBL_MAX;
  m->ef.upper_bound = +DBL_MAX;
}

// void liknorm_set_poisson(LikNormMachine *machine, double k) {
//   set_expfam(machine, POISSON, k, 1, -logfactorial(k),
//              "poisson");
// }
// void liknorm_set_exponential(LikNormMachine *machine, double x) {
//   set_expfam(machine, EXPONENTIAL, x, 1, 0, "exponential");
// }
// void liknorm_set_gamma(LikNormMachine *machine, double x, double a) {
//   set_expfam(machine, GAMMA, x, 1 / a, 0, "gamma");
// }
// void liknorm_set_geometric(LikNormMachine *machine, double x) {
//   set_expfam(machine, GEOMETRIC, x, 1, 0, "geometric");
// }

void liknorm_set_prior(LikNormMachine *machine, double tau, double eta) {
  const double tau_min = 2 * sqrt(DBL_EPSILON);
  tau = fmax(tau, tau_min);
  machine->normal.eta = eta;
  machine->normal.tau = tau;
  machine->normal.log_tau = log(tau);
}
