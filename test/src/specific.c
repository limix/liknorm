#include "hope.h"
#include "liknorm/liknorm.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#if (_MSC_VER <= 1500)
#include <float.h>
#define isnan(x) _isnan(x)
#define isfinite(x) _finite(x)
#endif
#endif

static void test_old_one(void);
static void test_marc(void);

int main()
{
    test_old_one();
    test_marc();
    return hope_status();
}

static void test_old_one(void)
{
    struct LikNormMachine *machine = liknorm_create_machine(500);
    double log_zeroth, mean, variance;

    double n = 142100;
    double k = 78155;
    liknorm_set_binomial(machine, k, n);
    liknorm_set_prior(machine, 1.108785906137200072407722473145,
                      8.231554676556697813794016838074);

    liknorm_integrate(machine, &log_zeroth, &mean, &variance);

    CLOSE(log_zeroth, -40.26009679144215169799);
    CLOSE(mean, 0.20089983974431713243);
    CLOSE2(variance, 0.00002843348172375942, 1e-08, 0);

    liknorm_destroy_machine(machine);
}

static void test_marc(void)
{
    struct LikNormMachine *machine = liknorm_create_machine(500);
    double log_zeroth, mean, variance;

    liknorm_set_bernoulli(machine, 0.5);
    double eta[] = {8.68315118780285e+19 * (0 / 100000.0),
                    8.68315118780285e+19 * (1 / 100000.0),
                    8.68315118780285e+19 * (337 / 100000.0),
                    8.68315118780285e+19 * (338 / 100000.0),
                    8.68315118780285e+19};

    double desired_log0[] = {-0.6969067887453335163883139102836139500141,
                             -13194139533312., 0., 0., 0.};
    double desired_mean[] = {0.0000000000000000989284694679729687520675,
                             26312579615912.1718750, 8867339309450631.000,
                             8893651729117382.00, 2631257911700760576.00};
    double desired_variance[] = {0.0300768633943054898571833888354376540519,
                                 DBL_EPSILON, DBL_EPSILON, DBL_EPSILON,
                                 DBL_EPSILON};

    for (unsigned i = 0; i < sizeof(desired_log0) / sizeof(double); ++i)
    {
        liknorm_set_prior(machine, 33.0, eta[i]);
        liknorm_integrate(machine, &log_zeroth, &mean, &variance);
        CLOSE(log_zeroth, desired_log0[i]);
        CLOSE(mean, desired_mean[i]);
        CLOSE(variance, desired_variance[i]);
    }

    liknorm_destroy_machine(machine);
}
