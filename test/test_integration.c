#include "liknorm/liknorm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#if (_MSC_VER <= 1500)
#define isnan(x) _isnan(x)
#define isfinite(x) _finite(x)
#endif
#endif

char *strdump(const char *str) {
  size_t i;
  char *rstr;
  for (i = 0;; ++i) {
    if (str[i] == '\0')
      break;
  }
  ++i;
  rstr = malloc(sizeof(char) * i);
  memcpy(rstr, str, i);
  return rstr;
}

typedef struct Line {
  double normal_mean;
  double normal_variance;
  char *likname;
  double y;
  double aphi;
  double mean;
  double variance;
  struct Line *next;
} Line;

Line *read_table() {
  FILE *stream = fopen("test/table.csv", "r");
  char line[65536];

  Line *l = 0;
  Line *root = 0;
  Line *last = root;
  char *token;
  char *str;

  if (stream == 0)
    return 0;

  str = fgets(line, 65536, stream);

  if (str == NULL)
    exit(1);

  while (fgets(line, 65536, stream)) {
    l = malloc(sizeof(Line));

    token = strtok(line, ",");
    l->normal_mean = atof(token);

    token = strtok(NULL, ",");
    l->normal_variance = atof(token);

    token = strtok(NULL, ",");
    l->likname = strdump(token);

    token = strtok(NULL, ",");
    l->y = atof(token);

    token = strtok(NULL, ",");
    l->aphi = atof(token);

    token = strtok(NULL, ",");
    l->mean = atof(token);

    token = strtok(NULL, ",");
    l->variance = atof(token);

    if (root == 0) {
      root = l;
      last = root;
    } else {
      last->next = l;
      last = l;
    }

    last->next = 0;
  }

  fclose(stream);

  return root;
}

int test_it(LikNormMachine *machine, Line *l, double *elapsed) {
  double log_zeroth, mean, variance;
  double eps = 1e-4;
  int ok;

  if (strcmp(l->likname, "bernoulli") == 0)
    liknorm_set_bernoulli(machine, l->y);

  if (strcmp(l->likname, "poisson") == 0)
    liknorm_set_poisson(machine, l->y);

  if (strcmp(l->likname, "binomial") == 0)
    liknorm_set_binomial(machine, l->y / l->aphi, 1 / l->aphi);

  if (strcmp(l->likname, "gamma") == 0)
    liknorm_set_gamma(machine, l->y, 1 / l->aphi);

  if (strcmp(l->likname, "exponential") == 0)
    liknorm_set_exponential(machine, l->y);

  if (strcmp(l->likname, "geometric") == 0)
    liknorm_set_geometric(machine, l->y);

  liknorm_set_prior(machine, 1 / l->normal_variance,
                    l->normal_mean / l->normal_variance);

  liknorm_integrate(machine, &log_zeroth, &mean, &variance);

  ok = fabs(mean - l->mean) < eps && fabs(variance - l->variance) < eps;
  ok = ok && isfinite(mean) && isfinite(variance);

  if (!ok) {
    printf("Test failed:\n");
    printf("name: %s\n", l->likname);
    printf("mean variance %g %g l->mean l->variance %g %g\n", mean, variance,
           l->mean, l->variance);
  }

  if (!ok)
    return 1;

  return 0;
}

int main() {
  Line *root = read_table();
  Line *l = root;
  int e;
  double elapsed = 0;

  LikNormMachine *machine = liknorm_create_machine(350);

  int i = 0;

  while (l != 0) {
    e = test_it(machine, l, &elapsed);

    if (e != 0)
      return 1;

    l = l->next;
    i++;
  }

  liknorm_destroy_machine(machine);

  return 0;
}
