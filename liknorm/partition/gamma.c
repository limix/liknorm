#include "gamma.h"

#include <assert.h>
#include <float.h>
#include <math.h>

double gamma_log_partition(const double theta) {
  assert(theta <= 0);
  return -log(fmax(DBL_EPSILON, -theta));
}

double gamma_log_partition_fderivative(const double theta) {
  assert(theta <= 0);
  return -log(fmax(DBL_EPSILON, -theta));
}

void gamma_log_partition_derivatives(const double theta, double *b0,
                                     double *logb1, double *logb2) {
  assert(theta <= 0);
  *b0 = -log(fmax(DBL_EPSILON, -theta));
  *logb1 = *b0;
  *logb2 = 2 * (*b0);
}
