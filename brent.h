#ifndef LIKNORM_BRENT_H
#define LIKNORM_BRENT_H

#include "func_base.h"

void liknorm_find_minimum(double *x0, double *fx0, liknorm_func_base *f,
                          void *args, double a, double b, double rtol,
                          double atol, int maxiter);

#endif
