#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "expfam.h"
#include "liknorm/liknorm.h"
#include "machine.h"
#include "normal.h"

void integrate_step(double si, double step, struct ExpFam *ef,
                    struct Normal *normal, double *log_zeroth, double *u,
                    double *v, double *A0, double *logA1, double *logA2,
                    double *diff);

void combine_steps(struct LikNormMachine *machine, double *log_zeroth,
                   double *mean, double *variance, double *left, double *right);

#endif
