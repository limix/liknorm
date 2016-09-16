#include "base.h"
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846


void integrate_step(double si, double step, double Ty,
                    log_partition *lp, Normal *normal, LogMoments *lm)
{
  double sii = si + step;
  double mi  = si / 2 + sii / 2;
  double eta = normal->eta;
  double tau = normal->tau;

  double A0, A1, A2;

  (*lp)(mi, &A0, &A1, &A2);

  double a = -A0 + A1 * mi - (A2 * mi * mi) / 2;
  double b = Ty + eta - A1 + A2 * mi;
  double c = -(tau + A2) / 2;

  *(lm->log_zeroth) = a - (b * b) / c + log(-PI / c) / 2;
}

void integrate(Interval *interval)
{}

int  main()
{
  Normal normal = { 0, 1 };

  printf("BLA BLA BLA");
}
