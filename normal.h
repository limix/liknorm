#ifndef LIKNORM_NORMAL_H
#define LIKNORM_NORMAL_H

struct Normal
{
    double eta;
    double log_tau;
    double tau;
};

/* Cumulative distribution function of the Normal distribution.
 */
double liknorm_cdf(const double x);

/* Log of the cumulative distribution function of the Normal distribution.
 */
double liknorm_logcdf(const double x);

/* Log of the probability distribution function of the Normal distribution.
 */
static inline double logpdf(const double x)
{
    return -(x * x) / 2 - 0.9189385332046726695409688545623794196;
}

#endif
