#include "liknorm/liknorm.h"
#include "benchmark.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* getfield(char *line, int num)
{
  line = strdup(line);
  char *tok;

  for (tok = strtok(line, ",");
       tok && *tok;
       tok = strtok(NULL, ",\n"))
  {
    if (!--num)
    {
      free(line);
      return strdup(tok);
    }
  }
  free(line);
  return NULL;
}

typedef struct Line
{
  double       normal_mean;
  double       normal_variance;
  char        *likname;
  double       y;
  double       aphi;
  double       mean;
  double       variance;
  struct Line *next;
} Line;

Line* read_table()
{
  FILE *stream = fopen("test/table.csv", "r");

  if (stream == 0) return 0;

  char line[65536];

  Line *l    = 0;
  Line *root = 0;
  Line *last = root;
  char *tmp;

  fgets(line, 65536, stream);

  while (fgets(line, 65536, stream))
  {
    l = malloc(sizeof(Line));

    tmp            = getfield(line, 1);
    l->normal_mean = atof(tmp);
    free(tmp);

    tmp                = getfield(line, 2);
    l->normal_variance = atof(tmp);
    free(tmp);

    l->likname = strdup(getfield(line, 3));

    tmp  = getfield(line, 4);
    l->y = atof(tmp);
    free(tmp);

    tmp     = getfield(line, 5);
    l->aphi = atof(tmp);
    free(tmp);

    tmp     = getfield(line, 6);
    l->mean = atof(tmp);
    free(tmp);

    tmp         = getfield(line, 7);
    l->variance = atof(tmp);
    free(tmp);

    if (root == 0)
    {
      root = l;
      last = root;
    } else {
      last->next = l;
      last       = l;
    }

    last->next = 0;
  }

  return root;
}

int test_it(LikNormMachine *machine, Line *l, double *elapsed)
{
  Normal normal;

  // printf("%s\n", l->likname);

  ExpFam ef =
  { l->y, l->aphi, log(l->aphi), get_log_partition(l->likname), 0, 0 };

  get_interval(l->likname, &(ef.left), &(ef.right));

  normal.tau     = 1 / l->normal_variance;
  normal.eta     = l->normal_mean / l->normal_variance;
  normal.log_tau = -log(l->normal_variance);

  double mean, variance;

  double start = get_time();
  integrate(machine, &ef, &normal, &mean, &variance);
  *elapsed += (get_time() - start);

  double eps = 1e-4;
  int    ok  = fabs(mean - l->mean) < eps && fabs(variance - l->variance) < eps;
  ok = ok && isfinite(mean) && isfinite(variance);

  if (!ok)
  {
    printf("%s %g %g %g %g %g %g\n",
           l->likname,
           l->normal_mean,
           l->normal_variance,
           ef.y,
           ef.aphi,
           mean - l->mean,
           variance - l->variance);
    printf("mean variance %g %g l->mean l->variance %g %g\n",
           mean,
           variance,
           l->mean,
           l->variance);
  }

  if (!ok) return 1;

  return 0;
}

int main()
{
  Line  *root = read_table();
  Line  *l    = root;
  int    e;
  double elapsed = 0;

  LikNormMachine *machine = create_liknorm_machine(350, 1e-7);

  int i = 0;

  while (l != 0)
  {
    e = test_it(machine, l, &elapsed);

    if (e != 0) return 1;

    l = l->next;
    i++;
  }

  destroy_liknorm_machine(machine);

  printf("Elapsed time: %.10f\n", elapsed / i);
  return 0;
}
