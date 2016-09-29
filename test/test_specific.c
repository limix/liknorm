#include "liknorm/liknorm.h"
#include "benchmark.h"
#include <math.h>
#include <stdio.h>


int main()
{
  LikNormMachine *machine = liknorm_create_machine(350);

  Normal normal;

  log_partition  *lp  = get_log_partition("binomial");
  log_partition0 *lp0 = get_log_partition0("binomial");

  ExpFam ef = { 0.93146417445482865, 1 / 321.0, log(1 / 321.0), lp, lp0, 0, 0,
                "binomial" };

  get_interval("binomial", &(ef.left), &(ef.right));

  normal.tau     = 4.1827910222757971;
  normal.eta     = -33.806632912457161;
  normal.log_tau = log(4.1827910222757971);

  double log_zeroth, mean, variance;

  liknorm_integrate(machine, &ef, &normal, &log_zeroth, &mean, &variance);

  if (fabs(log_zeroth + 1289.78232538217889668885618448) > 1e-7) return 1;

  if (fabs(mean + 4.66319474116641163874419362401) > 1e-7) return 1;

  if (fabs(variance - 1.25962109436272839957382529974e-05) > 1e-7) return 1;

  ef.y        = 0;
  ef.aphi     = 1. / 266;
  ef.log_aphi = log(ef.aphi);

  // normal.tau  = 3.97316;
  normal.tau = 1.0;

  // normal.eta = 1.44005e+07;

  const double total   = 100;
  const double inicial = 0;
  const double final   = 7000;

  for (int i = 0; i < total; ++i)
  {
    double eta = inicial + (final - inicial) * (i / total);
    normal.eta = eta;
    liknorm_integrate(machine, &ef, &normal, &log_zeroth, &mean, &variance);
    printf("%0.30g,%0.30g,%0.30g,%0.30g\n", normal.eta, log_zeroth, mean,
           variance);
  }

  //
  // printf("eta %g %.30g\n", normal.eta, log_zeroth);
  // printf("%.30g\n",        mean);
  // printf("%.30g\n",        variance);

  liknorm_destroy_machine(machine);

  return 0;
}
