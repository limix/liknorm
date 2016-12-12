#include "liknorm/liknorm.h"

int main() {
  LikNormMachine *machine = liknorm_create_machine(100000);

  double tau = 4.1827910222757971;
  double eta = -33.806632912457161;

  liknorm_set_prior(machine, tau, eta);

  double k = 3;
  double n = 10;
  liknorm_set_binomial(machine, k, n);

  double log_zeroth, mean, variance;
  liknorm_integrate(machine, &log_zeroth, &mean, &variance);

  // log_partition  *lp  = get_log_partition("binomial");
  // log_partition0 *lp0 = get_log_partition0("binomial");
  // log_partition0 *lp1 = get_log_partition1("binomial");

  // ExpFam ef =
  // { 0.93146417445482865,
  //   1 / 321.0,          log(
  //     1 / 321.0),       lp,
  //   lp0,
  //   lp1,
  //   0,
  //   0,
  //   "binomial" };

  // get_interval("binomial", &(ef.left), &(ef.right));

  // // if (fabs(log_zeroth + 1289.78232538217889668885618448) > 1e-7) return 1;
  // //
  // // if (fabs(mean + 4.66319474116641163874419362401) > 1e-7) return 1;
  // //
  // // if (fabs(variance - 1.25962109436272839957382529974e-05) > 1e-7) return
  // 1;
  //
  //
  // ef.y        = 0;
  // ef.aphi     = 1. / 266;
  // ef.log_aphi = log(ef.aphi);
  //
  // normal.tau = 1.0;
  //
  // // // normal.eta = 1.44005e+07;
  // //
  // // normal.tau = 1 / (10.0 * 10.0);
  // // normal.eta = 700 * normal.tau;
  // //
  // // liknorm_integrate(machine, &ef, &normal, &log_zeroth, &mean, &variance);
  // // return 0;
  //
  // //
  // // printf("%0.30g, %0.30g, %0.30g\n", log_zeroth, mean, variance);
  //
  //
  // const int total   = 100;
  // const int inicial = 0;
  // const int final   = 5000;
  //
  // // const double means[] = { 1, 100, 1000, 5000, 5001 };
  // // const double means[] = { 1, 100, 1000, 5000, 5001 };
  //
  // for (int i = 0; i < total; ++i)
  // {
  //   int eta = inicial + (final - inicial) * (((double)i) / total);
  //
  //   // double eta = means[i];
  //   normal.eta = eta;
  //   liknorm_integrate(machine, &ef, &normal, &log_zeroth, &mean, &variance);
  //   printf("%d,%0.30g,%0.30g,%0.30g\n", eta, log_zeroth, mean,
  //          variance);
  // }
  //
  // //
  // // printf("eta %g %.30g\n", normal.eta, log_zeroth);
  // // printf("%.30g\n",        mean);
  // // printf("%.30g\n",        variance);
  //
  liknorm_destroy_machine(machine);

  return 0;
}
