#ifndef GFUNC_H
#define GFUNC_H

#include "machine.h"

double g_function(double x, ExpFam *ef, Normal *normal);
double g_function_func_base(double x, void *args);
double g_derivative(double x, ExpFam *ef, Normal *normal);

#endif
