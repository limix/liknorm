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

int main() {
    struct LikNormMachine *machine = liknorm_create_machine(500);
    double log_zeroth, mean, variance;
    double eps = 1e-4;
    int ok;

    liknorm_set_binomial(machine, 784, 1421);

    liknorm_set_prior(machine, 1.108785906137200072407722473145,
                      8.231554676556697813794016838074);

    liknorm_integrate(machine, &log_zeroth, &mean, &variance);

    ok = fabs(mean - 0.78159506985534710211) < eps &&
         fabs(variance - 0.00002868685045742669) < eps;
    ok = ok && isfinite(mean) && isfinite(variance);

    liknorm_destroy_machine(machine);

    if (!ok)
        return 1;

    return 0;
}
