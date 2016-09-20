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

  LikNormMachine *machine = create_liknorm_machine(10, 1e-7);

  double N  = 10;
  double K  = 5;
  ExpFam ef = { K, &binomial_log_partition, &N };

  double mean, variance;

  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean) > 1e-10) return 1;

  destroy_liknorm_machine(machine);

  double mu;
  double var;

  machine    = create_liknorm_machine(1000, 1e-7);
  N          = 3;
  K          = 1;
  mu         = 0;
  var        = 1;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 0.301985252166403428386587393106) > 1e-7) return 1;

  N          = 9;
  K          = 2;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 0.838185833081513509412729945325) > 1e-7) return 1;

  N          = 9;
  K          = 8;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean - 1.214643858827865630090059312352) > 1e-7) return 1;

  N          = 9;
  K          = 8;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = 1.2;
  var        = 1.0;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean - 1.729607945703560689665323479858) > 1e-7) return 1;

  N          = 9;
  K          = 8;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = 1.2;
  var        = 0.3;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean - 1.442314118365103814412009342050) > 1e-7) return 1;

  N          = 29;
  K          = 2;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = 1.2;
  var        = 2.1;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 2.074258852713208423068635966047) > 1e-7) return 1;

  N          = 29;
  K          = 2;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = -9.2;
  var        = 2.1;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 5.482916813361223162814894749317) > 1e-7) return 1;

  N          = 5;
  K          = 0;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = -9.2;
  var        = 2.1;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 9.202995428625960983026743633673) > 1e-7) return 1;

  N          = 5;
  K          = 5;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = -1.2;
  var        = 2.1;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean - 1.332358223280064590809956825979) > 1e-7) return 1;


  N          = 5;
  K          = 0;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = -101.2;
  var        = 2.1;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 101.200000000255741383625718299299) > 1e-7) return 1;

  N          = 5;
  K          = 5;
  ef.Ty      = K;
  ef.lp      = &binomial_log_partition;
  ef.lp_data = &N;
  mu         = -101.2;
  var        = 2.1;
  normal.tau = 1 / var;
  normal.eta = mu / var;
  integrate(machine, ef, normal, &mean, &variance);

  if (fabs(mean + 90.713060785598017332631570752710) > 1e-7) return 1;

  destroy_liknorm_machine(machine);

  return 0;
}
