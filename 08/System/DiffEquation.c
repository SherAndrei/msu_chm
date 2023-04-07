#include "Error.h"
#include "Root.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// See Newton.pdf 2a

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
         "Usage: %s [-e <1e-10>] [-m <10>] [-X <1.0>] [-a <0>] [-b <0>]\n"
         "\n"
         "DESCRIPTION:\n"
         "\tFind numerical solution of differential equation y''=f(x,y)\n"
         "\tfor x in (0, X), y(0) = a, y(X) = b, found by\n"
         "\tsolving system of m non-linear algebraic equations\n"
         "\twith Newton's method.\n"
         "\n"
         "OPTIONS:\n"
         "\t-h            -- see help\n"
         "\t-e ( = 1e-10) -- precision of the solution, eps > 0\n"
         "\t-m ( = 10   ) -- amounf of equations in the system (unsigned)\n"
         "\t-X ( = 1.0  ) -- right bound for x, X > 0\n"
         "\t-a ( = 0    ) -- value of y(0)\n"
         "\t-b ( = 0    ) -- value of y(X)\n"
         "\t-v            -- verbose\n",
        argv0);
  return error;
}

static int ParseDouble(const char* from, double* value, const char* hint)
{
  if (sscanf(from, "%lf", value) != 1) {
    fprintf(stderr, "error: expected double for %s, got \"%s\"\n", hint, from);
    return InputError;
  }
  return Success;
}

static void InitialApproximation(double *x, unsigned m) {
  for (unsigned i = 0; i < m; i++) {
    x[i] = 0.;
  }
}

static double RightPart(double x, double y) { return cos(x) * sin(y); }

static double RightPartDerivativeByY(double x, double y) { return cos(x) * cos(y); }

static double X = 1.;
static double a = 0.;
static double b = 0.;

// amount of algebraic equations is N - 1
// amounf of equations required for Newton is m
// m == N - 1

static void F(double *y, const double *x, unsigned m) {
  const double h = X / (m + 1);
  const double divide_by_hh = 1. / h / h;

  y[0] = (x[1] - 2. * x[0] + a) * divide_by_hh - RightPart(h, x[0]);
  y[m - 1] = (b - 2. * x[m - 1] + x[m - 2]) * divide_by_hh - RightPart(X - h, x[m - 1]);
  for (unsigned k = 1; k < m - 1; k++) {
    y[k] = (x[k + 1] - 2. * x[k] + x[k - 1]) * divide_by_hh - RightPart((k + 1) * h, x[k]);
  }
}

#define J(i, j) jacobian[(i)*m + (j)]

static void dF(double *jacobian, const double *x, unsigned m) {
  const double h = X / (m + 1);
  const double main_diagonal = -2. / h / h;
  const double subdiagonal = 1. / h / h;
  const double supradiagonal = subdiagonal;

  for (unsigned i = 0; i < m * m; i++) {
    jacobian[i] = 0.;
  }

  J(0, 0) = main_diagonal - RightPartDerivativeByY(h, x[0]);
  J(0, 1) = supradiagonal;
  J(m - 1, m - 1) = main_diagonal - RightPartDerivativeByY(X - h, x[m - 1]);
  J(m - 1, m - 2) = subdiagonal;
  for (unsigned k = 1; k < m - 1; k++) {
    J(k, k - 1) = subdiagonal;
    J(k, k + 0) = main_diagonal - RightPartDerivativeByY((k + 1) * h, x[k]);
    J(k, k + 1) = supradiagonal;
  }
}

static inline void PrintSolution(const double *y, unsigned m) {
  const double h = X / (m + 1);
  printf("%20.14lf %20.14lf\n", 0., a);
  for (unsigned k = 0; k < m; k++) {
    printf("%20.14lf %20.14lf\n", (k + 1) * h, y[k]);
  }
  printf("%20.14lf %20.14lf\n", X, b);
}

int main(int argc, char * const *argv) {
  double eps = 1e-10;
  unsigned m = 10u;
  double *x;
  int root_search_result;
  int opt;
  int verbose = 0;

  while ((opt = getopt(argc, argv, "va:b:e:m:X:o:h")) != -1) {
    switch (opt) {
    case 'a':
      if (ParseDouble(optarg, &a, "a") != Success)
        return InputError;
      break;
    case 'b':
      if (ParseDouble(optarg, &b, "b") != Success)
        return InputError;
      break;
    case 'e':
      if (ParseDouble(optarg, &eps, "eps") != Success)
        return InputError;
      if (eps <= 0) {
        fprintf(stderr, "error: expected eps > 0, got \"%lf\"\n", eps);
        return InputError;
      }
      break;
    case 'X':
      if (ParseDouble(optarg, &X, "X") != Success)
        return InputError;
      if (X <= 0) {
        fprintf(stderr, "error: expected X > 0, got \"%lf\"\n", X);
        return InputError;
      }
      break;
    case 'm':
      if (sscanf(optarg, "%u", &m) != 1) {
        fprintf(stderr, "error: expected unsigned for m, got \"%s\"\n", optarg);
        return InputError;
      }
      if (m == 0) {
        fprintf(stderr, "error: expected m > 0\n");
        return InputError;
      }
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
      return Usage(argv[0], Success);
    default:
      return Usage(argv[0], InputError);
    }
  }

  x = malloc(m * sizeof(*x));
  if (!x) {
    fprintf(stderr, "error: not enough memory\n");
    return NotEnoughMemory;
  }

  InitialApproximation(x, m);
  if (verbose) {
    fprintf(stderr, "Initial approximation:\n");
    for (unsigned k = 0; k < m; k++) {
      fprintf(stderr, "%20.14lf\n", x[k]);
    }
    fputc('\n', stderr);
  }

  root_search_result = Root(x, F, dF, m, eps, verbose);

  if (root_search_result == 0) {
    PrintSolution(x, m);
  }

  free(x);
  return root_search_result;
}
