#include "liknorm.h"
#include <stdio.h>
#include <math.h>

double logaddexp(double x, double y)
{
  double tmp = x - y;

  if (x == y) return x + M_LN2;

  if (tmp > 0) return x + log1p(exp(-tmp));
  else if (tmp <= 0) return y + log1p(exp(tmp));

  return tmp;
}

// A(X) = N log(1 + e^x)
// A'(X) = N e^x / (1 + e^x)
// A''(X) = N e^x / (1 + e^x)^2
void binomial_log_partition(double  x,
                            double *A0,
                            double *logA1,
                            double *logA2)
{
  double N = 1;

  *A0    = N * logaddexp(0, x);
  *logA1 = log(N) + x - logaddexp(0, x);
  *logA2 = log(N) + x - 2 * logaddexp(0, x);
}

int main()
{
  Normal normal = { 0, 1 };
  double A0, logA1, logA2;

  binomial_log_partition(1.2, &A0, &logA1, &logA2);

  if (fabs(A0 - 1.46328246733803113422) > 1e-10) return 1;

  if (fabs(logA1 - -0.26328246733803117863) > 1e-10) return 1;

  if (fabs(logA2 - -1.72656493467606231285) > 1e-10) return 1;

  LikNormMachine *machine = create_liknorm_machine(10, 1e-7);

  // ExpFam ef               = { 5, &binomial_log_partition };
  // Normal normal           = { 0, 1 };
  //
  // integrate(machine, ef, normal);

  destroy_liknorm_machine(machine);

  return 0;
}
