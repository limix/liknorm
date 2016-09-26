#include "liknorm_python.h"
#include "liknorm.h"
#include <math.h>

static LikNormMachine *machine;

const char *liknames[] = { "binomial",    "bernoulli",    "poisson",    "gamma",
                           "exponential", "geometric" };


void py_create_liknorm_machine(int n, double precision)
{
  machine = create_liknorm_machine(350, 1e-7);
}

void py_destroy_liknorm_machine()
{
  destroy_liknorm_machine(machine);
}

void py_integrate(int likname_id, double *y, double *aphi,
                  double *normal_tau, double *normal_eta, int n,
                  double *mean, double *variance)
{
  Normal normal;
  ExpFam ef;
  log_partition *lp = get_log_partition(liknames[likname_id]);

  get_interval(liknames[likname_id], &(ef.left), &(ef.right));

  for (int i = 0; i < n; ++i)
  {
    ef.y           = y[i];
    ef.aphi        = aphi[i];
    ef.log_aphi    = log(ef.aphi);
    ef.lp          = lp;
    normal.tau     = normal_tau[i];
    normal.eta     = normal_eta[i];
    normal.log_tau = -log(normal.tau);

    integrate(machine, &ef, &normal, mean + i, variance + i);
  }
}
