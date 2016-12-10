#ifndef NORMAL_H
#define NORMAL_H

// #include <stddef.h>
// #include <float.h>
// #include <math.h>
//
// #include "constants.h"

typedef struct
{
  double eta;
  double log_tau;
  double tau;
} Normal;


/* Cumulative distribution function of the Normal distribution.
 */
double cdf(double x);

/* Log of the cumulative distribution function of the Normal distribution.
 */
double logcdf(double x);

/* Log of the probability distribution function of the Normal distribution.
 */
inline
double logpdf(double x)
{
  return - (x*x)/2 - 0.9189385332046726695409688545623794196;
}

#endif
