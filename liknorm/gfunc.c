#include "gfunc.h"
#include "expfam.h"
#include "normal.h"
#include <assert.h>
#include <math.h>

double g_function(double x, ExpFam *ef, Normal *normal)
{
  const double a = x * (ef->y / ef->aphi + normal->eta);
  const double b = (normal->tau * x * x) / 2;
  const double c = ef->lp(x) / ef->aphi;
  assert(isfinite(a));
  assert(isfinite(b));
  assert(isfinite(c));

  return (a - b) - c;
}

double g_function_func_base(double x, void *args)
{
  void **args_ = args;
  ExpFam *ef = args_[0];
  Normal *normal = args_[1];

  return g_function(x, ef, normal);
}

double g_derivative(double x, ExpFam *ef, Normal *normal)
{
  return ef->y / ef->aphi + normal->eta - x * normal->tau -
         exp(ef->lp1(x)) / ef->aphi;
}
