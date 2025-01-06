#ifndef LIKNORM_BRENT_H
#define LIKNORM_BRENT_H

#include "func_base.h"

void liknorm_find_minimum(double *x0, double *fx0, func_base *f, void *args,
                          double a, double b, double rtol, double atol,
                          int maxiter);

static inline double neg_func_base(double x, void *args)
{
    void **args_ = args;
    func_base *fb = (func_base *)args_[0];
    return -(*fb)(x, args_[1]);
}

static inline void find_maximum(double *x0, double *fx0, func_base *f,
                                void *args, double a, double b, double rtol,
                                double atol, int maxiter)
{
    void *args_[] = {(void *)f, args};
    liknorm_find_minimum(x0, fx0, &neg_func_base, args_, a, b, rtol, atol,
                         maxiter);
    *fx0 = -(*fx0);
}

#endif
