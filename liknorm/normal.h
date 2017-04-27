#ifndef NORMAL_H
#define NORMAL_H

typedef struct {
  double eta;
  double log_tau;
  double tau;
} Normal;

/* Cumulative distribution function of the Normal distribution.
 */
double cdf(const double x);

/* Log of the cumulative distribution function of the Normal distribution.
 */
double logcdf(const double x);

/* Log of the probability distribution function of the Normal distribution.
 */
static inline double logpdf(const double x) {
  return -(x * x) / 2 - 0.9189385332046726695409688545623794196;
}

#endif
