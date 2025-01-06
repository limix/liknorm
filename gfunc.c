#include "gfunc.h"
#include "expfam.h"
#include "normal.h"

double liknorm_g_function(double x, struct ExpFam *ef, struct Normal *normal)
{
    const double a = x * (ef->y / ef->a + normal->eta);
    const double b = (normal->tau * x * x) / 2;
    const double c = ef->lp(x) / ef->a;

    return (a - b) - c;
}

double liknorm_g_function_func_base(double x, void *args)
{
    void **args_ = args;
    struct ExpFam *ef = args_[0];
    struct Normal *normal = args_[1];

    return liknorm_g_function(x, ef, normal);
}
