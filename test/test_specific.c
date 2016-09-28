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

  ExpFam ef = { 0.93146417445482865, 1 / 321.0, log(1 / 321.0), lp, lp0, 0, 0 };

  get_interval("binomial", &(ef.left), &(ef.right));

  normal.tau     = 4.1827910222757971;
  normal.eta     = -33.806632912457161;
  normal.log_tau = log(4.1827910222757971);

  double log_zeroth, mean, variance;

  liknorm_integrate(machine, &ef, &normal, &log_zeroth, &mean, &variance);

  if (fabs(log_zeroth + 1289.78232538217889668885618448) > 1e-7) return 1;

  if (fabs(mean + 4.66319474116641163874419362401) > 1e-7) return 1;

  if (fabs(variance - 1.25962109436272839957382529974e-05) > 1e-7) return 1;


  liknorm_destroy_machine(machine);

  return 0;
}
