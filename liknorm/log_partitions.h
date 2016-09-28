#ifndef LOG_PARTITIONS_H
#define LOG_PARTITIONS_H

static inline double binomial_log_partition0(double theta)
{
  return theta + log1p(exp(-theta));
}
// \phi = N
// a(\phi) = 1/\phi
// b(\theta) = log(1 + e^\theta)
// b'(\theta) = e^\theta / (1 + e^\theta)
// b''(\theta) = e^\theta / (1 + e^\theta)^2
static inline
void binomial_log_partition(double  theta,
                            double *b0,
                            double *logb1,
                            double *logb2)
{
  *b0 = binomial_log_partition0(theta);
  *logb1 = theta - *b0;
  *logb2 = theta - 2 * (*b0);
}

static inline double bernoulli_log_partition0(double theta)
{
  return binomial_log_partition0(theta);
}
static inline
void bernoulli_log_partition(double  theta,
                             double *b0,
                             double *logb1,
                             double *logb2)
{
  binomial_log_partition(theta, b0, logb1, logb2);
}

static inline double poisson_log_partition0(double theta)
{
    return exp(theta);
}
static inline
void poisson_log_partition(double  theta,
                           double *b0,
                           double *logb1,
                           double *logb2)
{
  *b0 = poisson_log_partition0(theta);
  *logb2 = *logb1 = theta;
}

static inline double gamma_log_partition0(double theta)
{
    return -log(-theta);
}
static inline
void gamma_log_partition(double  theta,
                         double *b0,
                         double *logb1,
                         double *logb2)
{
  *b0 = gamma_log_partition0(theta);
  *logb1 = *b0;
  *logb2 = 2 * (*b0);
}

static inline double exponential_log_partition0(double theta)
{
    return gamma_log_partition0(theta);
}
static inline
void exponential_log_partition(double  theta,
                               double *b0,
                               double *logb1,
                               double *logb2)
{
  return gamma_log_partition(theta, b0, logb1, logb2);
}

static inline double geometric_log_partition0(double theta)
{
    return -log1p(-exp(theta));
}
static inline
void geometric_log_partition(double  theta,
                             double *b0,
                             double *logb1,
                             double *logb2)
{
    *b0 = geometric_log_partition0(theta);
    *logb1 = theta + *b0;
    *logb2 = theta + 2 * (*b0);
}

log_partition* get_log_partition(const char *name)
{
  if (strcmp(name, "binomial") == 0) return binomial_log_partition;

  if (strcmp(name, "bernoulli") == 0) return bernoulli_log_partition;

  if (strcmp(name, "poisson") == 0) return poisson_log_partition;

  if (strcmp(name, "gamma") == 0) return gamma_log_partition;

  if (strcmp(name, "exponential") == 0) return exponential_log_partition;

  if (strcmp(name, "geometric") == 0) return geometric_log_partition;

  return 0;
}

log_partition0* get_log_partition0(const char *name)
{
  if (strcmp(name, "binomial") == 0) return binomial_log_partition0;

  if (strcmp(name, "bernoulli") == 0) return bernoulli_log_partition0;

  if (strcmp(name, "poisson") == 0) return poisson_log_partition0;

  if (strcmp(name, "gamma") == 0) return gamma_log_partition0;

  if (strcmp(name, "exponential") == 0) return exponential_log_partition0;

  if (strcmp(name, "geometric") == 0) return geometric_log_partition0;

  return 0;
}

void get_interval(const char *name, double *left, double *right)
{
  *left  = -DBL_MAX;
  *right = +DBL_MAX;

  if (strcmp(name, "gamma") == 0) *right = -1e-15;

  if (strcmp(name, "exponential") == 0) *right = -1e-15;

  if (strcmp(name, "geometric") == 0) *right = -1e-15;
}

#endif /* end of include guard: LOG_PARTITIONS_H */
