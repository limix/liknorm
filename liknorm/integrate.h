#ifndef INTEGRATE_H
#define INTEGRATE_H

typedef struct ExpFam ExpFam;
typedef struct Normal Normal;

void integrate_step(double si, double step, ExpFam *ef, Normal *normal,
                    double *log_zeroth, double *u, double *v, double *A0,
                    double *logA1, double *logA2, double *diff);

#endif
