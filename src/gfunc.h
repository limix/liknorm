#ifndef GFUNC_H
#define GFUNC_H

#include "machine.h"

double g_function(double x, struct ExpFam *ef, struct Normal *normal);
double g_function_func_base(double x, void *args);

#endif
