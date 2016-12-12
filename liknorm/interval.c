#include "interval.h"
#include "expfam.h"
#include "normal.h"
#include "gfunc.h"
#include "optimizer/optimizer.h"
#include <math.h>
#include <assert.h>

static const double times_std = 7;
static const double reps = 1e-5;
static const double aeps = 1e-5;
static const double maxiter = 100;

static inline void find_first_interval(ExpFam *ef, Normal *normal, double *a,
                                       double *b) {
  double std = sqrt(1 / normal->tau);
  double mu = normal->eta / normal->tau;

  /* initial reasonable interval */
  *a = mu - times_std * std;
  *b = mu + times_std * std;

  *a = fmax(*a, ef->lower_bound);
  *b = fmin(*b, ef->upper_bound);

  const double smallest_step = fmax(fabs(*a), fabs(*b)) * reps + aeps;

  if (*b - *a < smallest_step) {
    if (*a - ef->lower_bound >= *b - ef->lower_bound)
      *a -= smallest_step;
    else
      *b += smallest_step;
  }
  assert(*b - *a >= smallest_step);
}

void find_interval(ExpFam *ef, Normal *normal, double *left, double *right) {

  double a, b;
  find_first_interval(ef, normal, &a, &b);
  void *args[] = {ef, normal};
  double fleft, fright;
  find_bracket(&g_function_func_base, args, a, b, ef->lower_bound,
               ef->upper_bound, left, right, &fleft, &fright);

  assert(*left < *right);
  assert(ef->lower_bound <= *left);
  assert(*right <= ef->upper_bound);
  assert(isfinite(*left) && isfinite(*right));
  assert(isfinite(fleft) && isfinite(fright));

  a = fmin(a, *left);
  b = fmax(b, *right);

  double xmax, fxmax;

  find_maximum(&xmax, &fxmax, &g_function_func_base, args, a, b, reps, aeps,
               maxiter);

  assert(isfinite(xmax));
  assert(isfinite(fxmax));

  assert(a <= xmax && xmax <= b);

  shrink_interval(ef, normal, &a, xmax, &b, fxmax);

  *left = a;
  *right = b;
}
