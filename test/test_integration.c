#include "liknorm/liknorm.h"
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
  if (strcmp(l->likname, "bernoulli") == 0)
    liknorm_set_bernoulli(machine, l->y);

  if (strcmp(l->likname, "poisson") == 0)
    liknorm_set_poisson(machine, l->y);

  if (strcmp(l->likname, "binomial") == 0)
    liknorm_set_binomial(machine, l->y / l->aphi, 1/l->aphi);

  if (strcmp(l->likname, "gamma") == 0)
    liknorm_set_gamma(machine, l->y, 1/l->aphi);

  if (strcmp(l->likname, "exponential") == 0)
    liknorm_set_exponential(machine, l->y);

  if (strcmp(l->likname, "geometric") == 0)
    liknorm_set_geometric(machine, l->y);

  liknorm_set_prior(machine, 1 / l->normal_variance,
                    l->normal_mean / l->normal_variance);

  double log_zeroth, mean, variance;
  liknorm_integrate(machine, &log_zeroth, &mean, &variance);

  double eps = 1e-4;

  int ok = fabs(mean - l->mean) < eps && fabs(variance - l->variance) < eps;
  ok = ok && isfinite(mean) && isfinite(variance);

  if (!ok)
  {
    printf("Test failed:\n");
    printf("name: %s\n", l->likname);
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

  LikNormMachine *machine = liknorm_create_machine(350);

  int i = 0;

  while (l != 0)
  {
    e = test_it(machine, l, &elapsed);

    if (e != 0) return 1;

    l = l->next;
    i++;
  }

  liknorm_destroy_machine(machine);

  return 0;
}
