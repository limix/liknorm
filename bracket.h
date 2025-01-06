#ifndef OPTIMIZER_BRACKET_H
#define OPTIMIZER_BRACKET_H

#include "func_base.h"

void liknorm_find_bracket(func_base *f, void *args, double a, double b,
                          double lower, double upper, double *left,
                          double *right, double *fleft, double *fright);

#endif
