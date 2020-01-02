#ifndef INTERVAL_C
#define INTERVAL_C

#include "interval.h"
#include "expfam.h"
#include "gfunc.h"
#include "hide.h"
#include "normal.h"
#include "optimizer/optimizer.h"
#include <math.h>

static const double times_std = 7;
static const double reps = 1e-5;
static const double aeps = 1e-5;
static const int maxiter = 100;

#ifndef DBL_TRUE_MIN
#define DBL_TRUE_MIN 4.9406564584124654E-324
#endif

static void find_first_interval(struct ExpFam *ef, struct Normal *normal, double *a,
                                double *b);
static double g_function_root(double x, void *args);
static void shrink_interval(struct ExpFam *ef, struct Normal *normal, double *a,
                            double xmax, double *b, double fxmax);

HIDE void find_interval(struct ExpFam *ef, struct Normal *normal, double *left,
                        double *right)
{

    double a, b;
    double xmax, fxmax;
    void *args[] = {ef, normal};
    double fleft, fright;

    find_first_interval(ef, normal, &a, &b);
    liknorm_find_bracket(&g_function_func_base, args, a, b, ef->lower_bound, ef->upper_bound,
                 left, right, &fleft, &fright);

    a = fmin(a, *left);
    b = fmax(b, *right);

    find_maximum(&xmax, &fxmax, &g_function_func_base, args, a, b, reps, aeps, maxiter);

    shrink_interval(ef, normal, &a, xmax, &b, fxmax);

    *left = a;
    *right = b;
}

static void find_first_interval(struct ExpFam *ef, struct Normal *normal, double *a,
                                double *b)
{
    double std = sqrt(1 / normal->tau);
    double mu = normal->eta / normal->tau;
    double smallest_step;

    /* initial reasonable interval */
    *a = mu - times_std * std;
    *b = mu + times_std * std;

    *a = fmax(*a, ef->lower_bound);
    *b = fmin(*b, ef->upper_bound);

    smallest_step = fmax(fabs(*a), fabs(*b)) * reps + aeps;

    if (*b - *a < smallest_step) {
        if (*a - ef->lower_bound >= *b - ef->lower_bound)
            *a -= smallest_step;
        else
            *b += smallest_step;
    }
}

static double g_function_root(double x, void *args)
{
    void **args_ = args;
    func_base *fb = (func_base *)args_[0];
    double *fxmax = args_[1];

    return *fxmax - (*fb)(x, args_[2]) + log(DBL_TRUE_MIN);
}

static void shrink_interval(struct ExpFam *ef, struct Normal *normal, double *a,
                            double xmax, double *b, double fxmax)
{
    void *args[] = {ef, normal};
    double fa = g_function_func_base(*a, args);
    double fb = g_function_func_base(*b, args);
    void *args_[] = {(void *)&g_function_func_base, &fxmax, args};

    if (fxmax - fa < log(DBL_TRUE_MIN)) {
        *a = zero(*a, xmax, 1e-5, &g_function_root, args_);
    }

    if (fxmax - fb < log(DBL_TRUE_MIN)) {
        *b = zero(*b, xmax, 1e-5, &g_function_root, args_);
    }
}

#endif
