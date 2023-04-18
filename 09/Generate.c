#include "Error.h"
#include "ExactSolution.h"

#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
         "Usage: %s N Lx1 Lx2\n"
         "DESCRIPTION:\n"
         "For a given rectangle [0, Lx1] x [0, Lx2] generate N random\n"
         "points using function from `ExactSolution.h`.\n"
         "\n"
         "OPTIONS:\n"
         "\tunsigned N - amount of points to generate, N > 1\n"
         "\tdouble Lx1 - upper bound of x1 coordinate\n"
         "\tdouble Lx2 - upper bound of x2 coordinate\n",
         argv0);
  return error;
}

static unsigned ParseToUnsigned(const char *param_name, char *from) {
  unsigned ret = 0;
  if (sscanf(from, "%u", &ret) != 1) {
    fprintf(stderr, "Error parsing %s has occured, terminate\n", param_name);
    exit(InputError);
  }
  from[0] = '\0';
  return ret;
}

static double ParseToDouble(const char *param_name, char *from) {
  double ret = 0.;
  if (sscanf(from, "%lf", &ret) != 1) {
    fprintf(stderr, "Error parsing %s has occured, terminate\n", param_name);
    exit(InputError);
  }
  from[0] = '\0';
  return ret;
}

int main(int argc, char *argv[]) {
  double Lx1 = 0.;
  double Lx2 = 0.;
  double x1, x2;
  unsigned N = 0;

  if (argc != 4)
    return Usage(argv[0], IncorrectUsage);

  N = ParseToUnsigned("N", argv[1]);
  if (N < 2) {
    fprintf(stderr, "error: N=%u is too small\n", N);
    return Usage(argv[0], InputError);
  }
  Lx1 = ParseToDouble("Lx1", argv[2]);
  Lx2 = ParseToDouble("Lx2", argv[3]);
  
  srand(42);
  printf("%u\n", N);
  for (unsigned i = 0u; i < N; i++) {
    x1 = (rand() / (RAND_MAX * 1.)) * Lx1;
    x2 = (rand() / (RAND_MAX * 1.)) * Lx2;
    printf("%20.15lf %20.15lf %20.15lf\n", x1, x2, ExactSolution(x1, x2));
  }

  return Success;
}
