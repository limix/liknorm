#ifndef GFUNC_H
#define GFUNC_H

struct ExpFam;
struct Normal;

double liknorm_g_function(double x, struct ExpFam *, struct Normal *);
double liknorm_g_function_func_base(double x, void *args);

#endif