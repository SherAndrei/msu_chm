#include "Error.h"
#include "Root.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
          "Usage: %s eps m [-v]\n"
          "DESCRIPTION:\n"
          "\tFind a numerical vector solution of system of `m` equations \n"
          "\tF(x)=0, x=(x^1,...,x^m)^T with an accuracy of `eps` using\n"
          "\tNewton's method.\n"
          "\n"
          "OPTIONS:\n"
          "\tdouble eps > 0 -- precision\n"
          "\tunsigned m > 0 -- amount of equations in system\n"
          "\t-v (optional) -- verbose\n",
          argv0);
  return error;
}

static inline void Print(FILE *out, const double *arr, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(out, "%20.14lf\n", arr[i]);
  }
}

#define J(i, j) jacobian[(i)*m + (j)]

static void InitialApproximation(double *x, unsigned m) {
  (void)m;
  x[0] = 0;
  x[1] = 0.5;
  x[2] = 1.;
}

static void F(double *f, const double *x, unsigned m) {
  (void)m;
  f[0] = pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2) - 10;
  f[1] = -pow(x[0], 2) + pow(x[1] - 3, 2) - x[2];
  f[2] = x[0] + 2. - x[1] / 3. + x[2];
}

static void dF(double *jacobian, const double *x, unsigned m) {
  (void)m;
  J(0, 0) = 2. * x[0];
  J(0, 1) = 2. * x[1];
  J(0, 2) = 2. * x[2];
  J(1, 0) = -2. * x[0];
  J(1, 1) = 2. * (x[1] - 3.);
  J(1, 2) = -1.;
  J(2, 0) = 1.;
  J(2, 1) = -1. / 3.;
  J(2, 2) = 1.;
}

int main(int argc, const char *argv[]) {
  double eps;
  unsigned m;
  double *x;
  int root_search_result;

  if (!(argc == 3 || argc == 4))
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &eps) == 1 && sscanf(argv[2], "%u", &m) == 1)) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return Usage(argv[0], InputError);
  }

  if (m == 0 || eps <= 0) {
    fprintf(stderr, "error: input should be greater than 0\n");
    return Usage(argv[0], InputError);
  }

  x = malloc(m * sizeof(*x));
  if (!x) {
    fprintf(stderr, "error: not enough memory\n");
    return NotEnoughMemory;
  }

  InitialApproximation(x, m);
  root_search_result = Root(x, F, dF, m, eps, argc == 4);

  if (root_search_result == 0) {
    Print(stdout, x, m);
  }

  free(x);
  return root_search_result;
}
