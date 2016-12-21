#include "liknorm/liknorm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
  LikNormMachine *machine = liknorm_create_machine(500);

  liknorm_set_binomial(machine, 784, 1421);

  liknorm_set_prior(machine, 1.108785906137200072407722473145,
                    8.231554676556697813794016838074);

  double log_zeroth, mean, variance;
  liknorm_integrate(machine, &log_zeroth, &mean, &variance);

  double eps = 1e-4;
  int ok = fabs(mean - 677.860089021189651248278096318245) < eps &&
           fabs(variance - 2.220446049250313080847263336182e-16) < eps;
  ok = ok && isfinite(mean) && isfinite(variance);

  liknorm_destroy_machine(machine);

  return 0;
}
