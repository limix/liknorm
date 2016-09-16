#ifndef MAIN_H
#define MAIN_H

typedef struct
{
  double left, right;
  int    n;
  double step;
} Interval;

typedef void log_partition (double  x,
                            double *A0,
                            double *A1,
                            double *A2);

typedef struct
{
  double *log_zeroth;
  double *log_first;
  double *log_second;
} LogMoments;

typedef struct
{
  double eta;
  double tau;
} Normal;

#endif /* ifndef MAIN_H */
