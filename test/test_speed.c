#include "liknorm.h"
#include "normal.h"
#include "logaddexp.h"
#include <stdio.h>
#include <math.h>

#ifdef WIN32

# include <windows.h>
double get_time()
{
  LARGE_INTEGER t, f;

  QueryPerformanceCounter(&t);
  QueryPerformanceFrequency(&f);
  return (double)t.QuadPart / (double)f.QuadPart;
}

#else /* ifdef WIN32 */

# include <sys/time.h>
# include <sys/resource.h>

double get_time()
{
  struct timeval  t;
  struct timezone tzp;

  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec * 1e-6;
}

#endif /* ifdef WIN32 */

int main()
{
  int nrep = 1000000;
  double r;
  double start = get_time();
  double tmp;

  for (int i = 0; i < nrep; i++)
  {
    r = exp(1.1 + 1. / i);

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for exp  : %.30f\n", (get_time() - start) / nrep);
  r += 1;

  start = get_time();

  for (int i = 0; i < nrep; i++)
  {
    r = log(1.1 + 1. / i);

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for log  : %.30f\n", (get_time() - start) / nrep);
  r += 1;

  start = get_time();

  for (int i = 0; i < nrep; i++)
  {
    r = logaddexp(-3.1 + 1. / i, 1.2 + 1. / (i + 1));

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for logae: %.30f\n", (get_time() - start) / nrep);
  r += 1;

  start = get_time();

  for (int i = 0; i < nrep; i++)
  {
    r = logcdf(-3.1 + 1. / i);

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for logcd: %.30f\n", (get_time() - start) / nrep);
  r += 1;

  start = get_time();

  for (int i = 0; i < nrep; i++)
  {
    r = logpdf(-3.1 + 1. / i);

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for logpd: %.30f\n", (get_time() - start) / nrep);
  r += 1;

  start = get_time();

  for (int i = 0; i < nrep; i++)
  {
    r = log1p(+3.1 + 1. / i);

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for log1p: %.30f\n", (get_time() - start) / nrep);
  r += 1;


  start = get_time();

  for (int i = 0; i < nrep; i++)
  {
    r = 1. / i + log1p(exp(-3.1 + 1. / i));

    if (r == 0.28723) printf("nao eh pra chegar aqui");
  }
  printf("Elapsed time for log1e: %.30f\n", (get_time() - start) / nrep);
  r += 1;

  return 0;
}
