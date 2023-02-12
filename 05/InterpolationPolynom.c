#include "Error.h"
#include "ExactSolution.h"
#include "GaussianElimination.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void Usage(const char *argv0) {
  printf("Usage: %s\n"
         "\tCalculate interpolation polynom by using input data.\n"
         "\n"
         "\tExpected input format:\n"
         "\t\tN\n"
         "\t\tx1\ty1\n"
         "\t\t...\t...\n"
         "\t\txN\tyN\n"
         "\n"
         "\tOutput format (3*N-2 rows):\n"
         "\t\tx1 \ty1 \tpolynom11\tdelta11\n"
         "\t\tx12\ty12\tpolynom12\tdelta12\n"
         "\t\tx13\ty13\tpolynom13\tdelta13\n"
         "\t\tx2 \ty2 \tpolynom21\tdelta21\n"
         "\t\tx22\ty22\tpolynom22\tdelta22\n"
         "\t\tx23\ty23\tpolynom23\tdelta23\n"
         "\t\t...\t...\t...\t...\n"
         "\t\txN \tyN \tpolynomN \tdeltaN \n",
         argv0);
}

static void FillVandermondeMatrix(const double *x, unsigned N, double *A) {
  for (unsigned i = 0u; i < N; i++)
    for (unsigned j = 0u; j < N; j++)
      A[i * N + j] = pow(x[i], j);
}

static int SolveWithGaussianElimination(const double *x, const double *y,
                                        unsigned N, double *a) {
  double *A = (double *)malloc(N * N * sizeof(double));
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

static double InterpolationResult(const double *a, double x, unsigned N) {
  double res = 0.;
  for (unsigned i = 0; i < N; i++)
    res += a[i] * pow(x, i);
  return res;
}

static void PrintSingleEntry(double x, const double *a, unsigned N) {
  const double exact = ExactSolution(x);
  const double result = InterpolationResult(a, x, N);
  fprintf(stdout, "%e %e %e %e\n", x, exact, result, fabs(result - exact));
}

static void PrintResult(const double *x, const double *a, unsigned N) {
  double step = 0.;
  for (unsigned i = 0; i + 1 < N; i++) {
    step = (x[i + 1] - x[i]) / 3;
    PrintSingleEntry(x[i] + 0 * step, a, N);
    PrintSingleEntry(x[i] + 1 * step, a, N);
    PrintSingleEntry(x[i] + 2 * step, a, N);
  }
  PrintSingleEntry(x[N - 1], a, N);
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  double *x = NULL;
  double *y = NULL;
  double *a = NULL;

  if (argc > 2) {
    Usage(argv[0]);
    return IncorrectUsage;
  }

  if (fscanf(stdin, "%u", &N) != 1) {
    fprintf(stderr, "Error parsing N from input\n");
    return InputError;
  }

  x = (double *)malloc(N * sizeof(double));
  y = (double *)malloc(N * sizeof(double));
  a = (double *)malloc(N * sizeof(double));
  if (!x || !y || !a) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    free(a);
    return NotEnoughMemory;
  }

  for (unsigned i = 0u; i < N; i++) {
    if (fscanf(stdin, "%lf%lf", x + i, y + i) != 2) {
      fprintf(stderr, "Error parsing input data\n");
      free(x);
      free(y);
      free(a);
      return InputError;
    }
  }

  SolveWithGaussianElimination(x, y, N, a);

  PrintResult(x, a, N);

  free(x);
  free(y);
  free(a);
  return Success;
}
