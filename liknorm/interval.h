#ifndef INTERVAL_H
#define INTERVAL_H

typedef struct ExpFam ExpFam;
typedef struct Normal Normal;

void find_interval(ExpFam *ef, Normal *normal, double *left, double *right);

#endif
