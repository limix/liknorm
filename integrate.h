#ifndef LIKNORM_INTEGRATE_H
#define LIKNORM_INTEGRATE_H

#include "expfam.h"
#include "liknorm.h"
#include "machine.h"
#include "normal.h"

void liknorm_integrate_step(double si, double step, struct ExpFam *,
                            struct Normal *, double *log_zeroth, double *u,
                            double *v, double *A0, double *logA1, double *logA2,
                            double *diff);

void liknorm_combine_steps(struct LikNormMachine *, double *log_zeroth,
                           double *mean, double *variance, double *left,
                           double *right);

#endif
