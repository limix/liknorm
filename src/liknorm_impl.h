#ifndef LIKNORM_IMPL_H
#define LIKNORM_IMPL_H

typedef struct
{
  double left, right;
  int    n;
  double step;
} Interval;

typedef struct
{
  double* log_zeroth;
  double* first;
  double* second;
} LogMoments;

#endif /* ifndef LIKNORM_IMPL_H */
