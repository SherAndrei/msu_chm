#include "Error.h"
#include "ExactSolution.h"
#include "GaussianElimination.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  printf("Usage: %s\n"
         "\tCalculate interpolation polynom of degree N-1 by using input data "
         "which consists of N points.\n"
         "\n"
         "\tExpected input format:\n"
         "\t\tN\n"
         "\t\tx1\ty1\n"
         "\t\t...\t...\n"
         "\t\txN\tyN\n"
         "\n"
         "\tOutput format (3*N-2 rows):\n"
         "\t\tx1 \ty1 \tcanonical11\tdelta11\tlagrangian11\tdelta11\n"
         "\t\tx12\ty12\tcanonical12\tdelta12\tlagrangian12\tdelta12\n"
         "\t\tx13\ty13\tcanonical13\tdelta13\tlagrangian13\tdelta13\n"
         "\t\tx2 \ty2 \tcanonical21\tdelta21\tlagrangian21\tdelta21\n"
         "\t\tx22\ty22\tcanonical22\tdelta22\tlagrangian22\tdelta22\n"
         "\t\tx23\ty23\tcanonical23\tdelta23\tlagrangian23\tdelta23\n"
         "\t\t...\t...\t...\t...\n"
         "\t\txN \tyN \tcanonicalN \tdeltaN \tlagrangianN \tdeltaN \n",
         argv0);
  return error;
}

static void FillVandermondeMatrix(const double *x, unsigned N, double *A) {
  for (unsigned i = 0u; i < N; i++)
    for (unsigned j = 0u; j < N; j++)
      A[i * N + j] = pow(x[i], j);
}

static int FindCanonicalCoefficients(const double *x, const double *y,
                                     unsigned N, double *a) {
  double *A = malloc(N * N * sizeof(*A));
  if (!A) {
    fprintf(stderr, "Not enough memory\n");
    return 4;
  }

  FillVandermondeMatrix(x, N, A);

  if (GaussMaxCol(A, N, y, a) < 0) {
    fprintf(stderr, "Matrix is degenerate!\n");
    free(A);
    return 1;
  }

  free(A);
  return 0;
}

static double CanonicalForm(const double *a, double x, unsigned N) {
  // Horner's method
  double res = a[N - 1] * x;
  for (unsigned i = N - 2; i > 0; --i) {
    res += a[i];
    res *= x;
  }
  return res + a[0];
}

static double PhiNumerator(unsigned i, unsigned N, double xi, const double *x) {
  double numerator = 1.;

  for (unsigned j = 0u; j < N; j++) {
    if (j == i)
      continue;

    numerator *= (xi - x[j]);
  }
  return numerator;
}

static double PhiDenominator(unsigned i, unsigned N, const double *x) {
  double denominator = 1.;

  for (unsigned j = 0u; j < N; j++) {
    if (j == i)
      continue;

    denominator *= (x[i] - x[j]);
  }
  return 1. / denominator;
}

static void FindLagrangianCoefficients(const double *x, const double *y,
                                       unsigned N, double *coefs) {
  for (unsigned i = 0; i < N; i++) {
    coefs[i] = y[i] * PhiDenominator(i, N, x);
  }
}

static double LagrangianForm(const double *x, const double *coefs, double xi,
                             unsigned N) {
  double res = 0.;
  for (unsigned i = 0; i < N; i++)
    res += coefs[i] * PhiNumerator(i, N, xi, x);
  return res;
}

static void PrintSingleEntry(double xi, const double *x,
                             const double *canonical_coefs,
                             const double *lagrangian_coefs, unsigned N) {
  const double exact = ExactSolution(xi);
  const double canonical = CanonicalForm(canonical_coefs, xi, N);
  const double lagrange = LagrangianForm(x, lagrangian_coefs, xi, N);
  fprintf(stdout, "%e %e %e %e %e %e\n", xi, exact, canonical,
          fabs(canonical - exact), lagrange, fabs(lagrange - exact));
}

static void PrintResult(const double *x, const double *canonical_coefs,
                        const double *lagrangian_coefs, unsigned N) {
  double step = 0.;
  for (unsigned i = 0; i + 1 < N; i++) {
    step = (x[i + 1] - x[i]) / 3;
    PrintSingleEntry(x[i] + 0 * step, x, canonical_coefs, lagrangian_coefs, N);
    PrintSingleEntry(x[i] + 1 * step, x, canonical_coefs, lagrangian_coefs, N);
    PrintSingleEntry(x[i] + 2 * step, x, canonical_coefs, lagrangian_coefs, N);
  }
  PrintSingleEntry(x[N - 1], x, canonical_coefs, lagrangian_coefs, N);
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  double *x = NULL;
  double *y = NULL;
  double *canonical_coefs = NULL;
  double *lagrangian_coefs = NULL;

  if (argc != 1)
    return Usage(argv[0], IncorrectUsage);

  if (fscanf(stdin, "%u", &N) != 1) {
    fprintf(stderr, "error: parsing N from input\n");
    return Usage(argv[0], InputError);
  }

  x = malloc(N * sizeof(*x));
  y = malloc(N * sizeof(*y));
  if (!x || !y) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    return NotEnoughMemory;
  }

  for (unsigned i = 0u; i < N; i++) {
    if (fscanf(stdin, "%lf%lf", x + i, y + i) != 2) {
      fprintf(stderr, "error: parsing input data\n");
      free(x);
      free(y);
      return Usage(argv[0], InputError);
    }
  }

  canonical_coefs = malloc(N * sizeof(*canonical_coefs));
  lagrangian_coefs = malloc(N * sizeof(*lagrangian_coefs));
  if (!canonical_coefs || !lagrangian_coefs) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    free(canonical_coefs);
    free(lagrangian_coefs);
    return NotEnoughMemory;
  }

  FindCanonicalCoefficients(x, y, N, canonical_coefs);
  FindLagrangianCoefficients(x, y, N, lagrangian_coefs);

  PrintResult(x, canonical_coefs, lagrangian_coefs, N);

  free(x);
  free(y);
  free(canonical_coefs);
  free(lagrangian_coefs);
  return Success;
}
