#ifndef LIKNORM_INTERVAL_H
#define LIKNORM_INTERVAL_H

#include "expfam.h"
#include "normal.h"

void liknorm_find_interval(struct ExpFam *ef, struct Normal *normal,
                           double *left, double *right);

#endif
