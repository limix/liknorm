#ifndef INTERVAL_H
#define INTERVAL_H

#include "expfam.h"
#include "normal.h"

void find_interval(struct ExpFam *ef, struct Normal *normal, double *left,
                   double *right);

#endif
