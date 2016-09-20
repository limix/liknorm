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
                            void   *lp_data,
                            double *A0,
                            double *logA1,
                            double *logA2,
                            double *sign)
{
  double N  = *((double *)lp_data);
  double ax = logaddexp(0, x);

  *A0    = N * ax;
  *logA1 = log(N) + x - ax;
  *logA2 = log(N) + x - 2 * ax;
  *sign  = +1;
}

int main()
{
  Normal normal = { 0, 1 };
  double A0, logA1, logA2, sign;

  // binomial_log_partition(1.2, &A0, &logA1, &logA2, &sign);

  // if (fabs(A0 - 1.46328246733803113422) > 1e-10) return 1;
  //
  // if (fabs(logA1 - -0.26328246733803117863) > 1e-10) return 1;
  //
  // if (fabs(logA2 - -1.72656493467606231285) > 1e-10) return 1;

  LikNormMachine *machine = create_liknorm_machine(10, 1e-7);

  double N  = 10;
  double K  = 5;
  ExpFam ef = { K, &binomial_log_partition, &N };

  double mean, variance;

  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean) > 1e-10) return 1;

  destroy_liknorm_machine(machine);

  N = 5;
  K = 2;

  // double mu  = 0.1;
  // double var = 1.2;

  // double lmom0   = -1.87923983187;
  // double mu_res  = -0.189211912705;
  // double var_res = 0.256669390778;

  return 0;
}
