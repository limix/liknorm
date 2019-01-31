#include "cass.h"
#include "liknorm.h"
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

int main()
{
    struct LikNormMachine *machine = liknorm_create_machine(500);
    double log_zeroth, mean, variance;

    double n = 142100;
    double k = 78155;
    liknorm_set_binomial(machine, k, n);
    liknorm_set_prior(machine, 1.108785906137200072407722473145,
                      8.231554676556697813794016838074);

    liknorm_integrate(machine, &log_zeroth, &mean, &variance);

    printf("%.20f\n%.20f\n%.20f\n", log_zeroth, mean, variance);

    assert_almost_equal(log_zeroth, -40.26009679144215169799);
    assert_almost_equal(mean, 0.20089983974431713243);
    assert_almost_equal(variance, 0.00002843348172375942);

    liknorm_destroy_machine(machine);

    return 0;
}
