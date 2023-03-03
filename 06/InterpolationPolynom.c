#include "Error.h"
#include "ExactSolution.h"
#include "GaussianElimination.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// from task:
// * input data consists n+2 points
// * polynom of degree n so there is n+1 coeff

// by me:
// * input data constists of N points
// * polynom of degree N-2 so there is N-1 coeff

static int Usage(const char *argv0, int error) {
  printf("Usage: %s\n"
         "\tCalculate interpolation polynom of degree N-2 by using input data which consists N points.\n"
         "\n"
         "\tExpected input format:\n"
         "\t\tN\n"
         "\t\tx1\ty1\n"
         "\t\t...\t...\n"
         "\t\txN\tyN\n"
         "\n"
         "\tOutput format (N rows):\n"
         "\t\tx1 \ty1 \tcanonical1\tdelta1\n"
         "\t\tx2 \ty2 \tcanonical2\tdelta2\n"
         "\t\t...\t...\t...\t...\n"
         "\t\txN \tyN \tcanonicalN \tdeltaN\n",
         argv0);
  return error;
}

static void FillVandermondeMatrixAndColumnWithH(const double *x, unsigned N, double *A) {
  for (unsigned i = 0u; i < N; i++) {
    for (unsigned j = 0u; j < N-1; j++) {
      A[i * N + j] = pow(x[i], j);
    }
    A[i * N + N - 1] = pow(-1., i);
  }
}

static int FindCanonicalCoefficients(const double *x, const double *y,
                                     unsigned N, double *a) {
  // last column filled up with (-1)^i
  double *A = (double *)malloc(N * N * sizeof(double));
  if (!A) {
    fprintf(stderr, "Not enough memory\n");
    return 4;
  }

  FillVandermondeMatrixAndColumnWithH(x, N, A);

  if (GaussMaxCol(A, N, y, a) < 0) {
    fprintf(stderr, "Matrix is degenerate!\n");
    free(A);
    return 1;
  }

  free(A);
  return 0;
}

static double CanonicalForm(const double *a, double x, unsigned N) {
  double res = 0.;
  for (unsigned i = 0; i < N; i++)
    res += a[i] * pow(x, i);
  return res;
}

static void PrintSingleEntry(double xi,
                             const double *canonical_coefs,
                             unsigned N) {
  const double exact = ExactSolution(xi);
  const double canonical = CanonicalForm(canonical_coefs, xi, N);
  fprintf(stdout, "%20e %20e %20e %20e\n", xi, exact, canonical,
          fabs(canonical - exact));
}

static void PrintResult(const double *x, const double *canonical_coefs_with_h,
                        unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    PrintSingleEntry(x[i], canonical_coefs_with_h, N);
  }
  fprintf(stdout, "\nh = %e\n", canonical_coefs_with_h[N - 1]);
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  double *x = NULL;
  double *y = NULL;
  double *canonical_coefs_with_h = NULL;

  if (argc != 1)
    return Usage(argv[0], IncorrectUsage);

  if (fscanf(stdin, "%u", &N) != 1) {
    fprintf(stderr, "error: parsing N from input\n");
    return Usage(argv[0], InputError);
  }

  x = (double *)malloc(N * sizeof(double));
  y = (double *)malloc(N * sizeof(double));
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
  // N-1 coeffs and h at the end
  canonical_coefs_with_h = (double *)malloc(N * sizeof(double));
  if (!canonical_coefs_with_h) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    free(canonical_coefs_with_h);
    return NotEnoughMemory;
  }

  FindCanonicalCoefficients(x, y, N, canonical_coefs_with_h);

  PrintResult(x, canonical_coefs_with_h, N);

  free(x);
  free(y);
  free(canonical_coefs_with_h);
  return Success;
}
