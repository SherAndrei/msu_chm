#include "Error.h"
#include "ExactSolution.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static void Usage(const char *argv0) {
  printf("Usage: %s N left right [OPTION]\n"
         "\tunsigned N - amount of points to generate, N > 1\n"
         "\tdouble left - left bound of the segment\n"
         "\tdouble right - right bound of the segment\n"
         "Options:\n"
         "\t-d\tgenerate x with distributed equally nodes\n"
         "\t-c\tgenerate x with Chebyshev nodes\n",
         argv0);
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

void GenerateChebyshevNodes(unsigned N, double left, double right, double *x) {
  const double add = 0.5 * (left + right);
  const double coef = 0.5 * (right - left);
  for (unsigned k = 0; k < N; k++)
    x[k] = add + coef * cos(M_PI * (2 * (k + 1) - 1) / (2. * N));
}

static void GenerateEquallyDistributedNodes(unsigned N, double left,
                                            double right, double *x) {
  double step = (right - left) / (N - 1.);
  x[0] = left;
  for (unsigned i = 0; i + 1 < N; i++)
    x[i + 1] = x[i] + step;
}

int main(int argc, char *argv[]) {
  double left_bound = 0.;
  double right_bound = 0.;
  double *x = NULL;
  unsigned N = 0;

  if (argc != 5) {
    Usage(argv[0]);
    return IncorrectUsage;
  }

  N = ParseToUnsigned("N", argv[1]);
  if (N < 2) {
    fprintf(stderr, "error: N=%u is too small\n", N);
    return IncorrectUsage;
  }
  left_bound = ParseToDouble("left bound", argv[2]);
  right_bound = ParseToDouble("right bound", argv[3]);

  x = (double *)malloc(sizeof(double) * N);
  if (!x) {
    fprintf(stderr, "Not enough memory\n");
    return NotEnoughMemory;
  }

  switch (getopt(argc, argv, "cd")) {
  case 'c':
    GenerateChebyshevNodes(N, left_bound, right_bound, x);
    break;
  case 'd':
    GenerateEquallyDistributedNodes(N, left_bound, right_bound, x);
    break;
  default:
    Usage(argv[0]);
    free(x);
    return InputError;
  }

  printf("%u\n", N);
  for (unsigned i = 0u; i < N; i++)
    printf("%e %e\n", x[i], ExactSolution(x[i]));

  free(x);
  return Success;
}
