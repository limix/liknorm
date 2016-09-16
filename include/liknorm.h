#ifndef LIKNORM_H
#define LIKNORM_H

typedef void log_partition (double  x,
                            double *A0,
                            double *A1,
                            double *A2);

typedef struct
{
  double eta;
  double tau;
} Normal;

#endif /* ifndef LIKNORM_H */
