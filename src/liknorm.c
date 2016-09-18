#include "liknorm_impl.h"
#include "liknorm.h"
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
  double mi  = si / 2 + sii / 2;
  double eta = normal.eta;
  double tau = normal.tau;

  double A0, A1, A2;

  (*ef.lp)(mi, &A0, &A1, &A2);

  double a = -A0 + A1 * mi - (A2 * mi * mi) / 2;
  double b = ef.Ty + eta - A1 + A2 * mi;
  double c = -(tau + A2) / 2;

  *(lm.log_zeroth) = a - (b * b) / c + log(-PI / c) / 2;
}

void integrate(LikNormMachine *machine, ExpFam ef, Normal normal)
{
  double left  = normal.eta / normal.tau - 10 * sqrt(1 / normal.tau);
  double right = normal.eta / normal.tau + 10 * sqrt(1 / normal.tau);

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
    lm.log_first  = machine->log_first + i;
    lm.log_second = machine->log_second + i;

    integrate_step(si, interval.step, ef, normal, lm);
  }
}

LikNormMachine* create_liknorm_machine(int n, double precision)
{
  LikNormMachine *machine = malloc(sizeof(LikNormMachine));

  machine->n          = n;
  machine->log_zeroth = malloc(n * sizeof(double));
  machine->log_first  = malloc(n * sizeof(double));
  machine->log_second = malloc(n * sizeof(double));
  machine->precision  = precision;

  return machine;
}

void destroy_liknorm_machine(LikNormMachine *machine)
{
  free(machine->log_zeroth);
  free(machine->log_first);
  free(machine->log_second);
  free(machine);
}
