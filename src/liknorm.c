#include "liknorm_impl.h"
#include "liknorm.h"
#include "logaddexp.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846


void integrate_step(double     si,
                    double     step,
                    ExpFam     ef,
                    Normal     normal,
                    LogMoments lm)
{
  double sii = si + step;

  printf("si sii %.10f %.10f\n", si, sii);
  double mi  = si / 2 + sii / 2;
  double eta = normal.eta;
  double tau = normal.tau;

  double A0, logA1, logA2, sign;
  (*ef.lp)(mi, &A0, &logA1, &logA2, &sign);
  printf("eta tau %.10f %.10f\n",                        eta, tau);
  printf("A0 logA1 logA2 sign %.10f %.10f %.10f %.1f\n", A0,  logA1, logA2,
         sign);

  double a = -A0 + exp(logaddexps(logA1, logA2, sign * mi, -(mi * mi) / 2));
  double bsign;
  double tmp = logaddexpss(logA1, logA2, -sign, mi, &bsign);
  double b   = ef.Ty + eta + copysign(exp(tmp), bsign);
  double c   = -(tau + exp(logA2)) / 2;
  double bc  = b / (2 * c);

  printf("a b c %.10f %.10f %.10f\n", a, b, c);

  *(lm.log_zeroth) = a - b * bc / 2 + log(PI) / 2 - log(-c) / 2;
  *(lm.first)      = -bc;
  printf("  lm.log_zeroth %.10f\n", *(lm.log_zeroth));
  printf("  lm.first      %.10f\n", *(lm.first));
  *(lm.second) = bc * bc - 1 / (2 * c);
}

void combine_steps(LikNormMachine *machine, double *mean, double *variance)
{
  (*mean)     = 0;
  (*variance) = 0;

  double total = logaddexp(machine->log_zeroth[0], machine->log_zeroth[1]);
  int    i;

  for (i = 2; i < machine->n; i++) total = logaddexp(total,
                                                     machine->log_zeroth[i]);

  for (i = 0; i < machine->n; i++)
  {
    (*mean) += machine->first[i] * exp(machine->log_zeroth[i] - total);
  }
}

void integrate(LikNormMachine *machine,
               ExpFam          ef,
               Normal          normal,
               double         *mean,
               double         *variance)
{
  double left  = normal.eta / normal.tau - 10 * sqrt(1 / normal.tau);
  double right = normal.eta / normal.tau + 10 * sqrt(1 / normal.tau);

  printf("left right %.10f %.10f\n", left, right);

  Interval   interval;
  LogMoments lm;

  interval.left  = left;
  interval.right = right;
  interval.n     = machine->n;
  interval.step  = (right - left) / interval.n;

  double si;

  for (int i = 0; i < interval.n; i++)
  {
    si            = interval.left + interval.step * i;
    lm.log_zeroth = machine->log_zeroth + i;
    lm.first      = machine->first + i;
    lm.second     = machine->second + i;

    integrate_step(si, interval.step, ef, normal, lm);
  }

  combine_steps(machine, mean, variance);

  printf("mean     %.10f\n", *mean);

  // printf("variance %.10f\n", *variance);
}

LikNormMachine* create_liknorm_machine(int n, double precision)
{
  LikNormMachine *machine = malloc(sizeof(LikNormMachine));

  machine->n          = n;
  machine->log_zeroth = malloc(n * sizeof(double));
  machine->first      = malloc(n * sizeof(double));
  machine->second     = malloc(n * sizeof(double));
  machine->precision  = precision;

  return machine;
}

void destroy_liknorm_machine(LikNormMachine *machine)
{
  free(machine->log_zeroth);
  free(machine->first);
  free(machine->second);
  free(machine);
}
