#include "Error.h"
#include "ExactSolution.h"
#include "GaussianElimination.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  printf("Usage: %s m\n"
         "DESCRIPTION:\n"
         "\tCalculate interpolation polynom of degree m < (N-1) by using input data which consists N > 1 points.\n"
	 "OPTIONS:\n"
	 "\tunsigned m -- degree of desired interpolation polynom\n"
         "\n"
         "INPUT FORMAT:\n"
         "\t\tN\n"
         "\t\tx1\ty1\n"
         "\t\t...\t...\n"
         "\t\txN\tyN\n"
         "\n"
         "OUTPUT FORMAT:\n"
	 "\tN rows, each constists starting nodes, exact solution, result of the interpolation polynom and their difference\n"
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

static int FindCanonicalCoefficients(const double *x, const double *basis,
                                     unsigned m, double *a) {
  // last column filled up with (-1)^i
  double *A = (double *)malloc(m * m * sizeof(double));
  if (!A) {
    fprintf(stderr, "Not enough memory for matrix\n");
    return 4;
  }

  FillVandermondeMatrixAndColumnWithH(x, m, A);

  if (GaussMaxCol(A, m, basis, a) < 0) {
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

// from N points select m+2 equally distributed points for basis
static void SelectBasis(const double* y, unsigned N, unsigned m, double* basis)
{
  unsigned step = N / ((m + 2u) - 1u);
  for (unsigned i = 0u, basis_i = 0u; i < N && basis_i < m; i += step, ++basis_i)
    basis[basis_i] = y[i];
}

// if maximum deviation is less than h, then H and H_pos is unchanged
static void FindMaximumDeviationAndItsPosition(const double* y, const double* canonical_coefs_with_h, unsigned m, double* H, unsigned* H_pos)
{
  const double h = canonical_coefs_with_h[m + 1];
  double current = 0.;
  for (unsigned i = 0; i < m + 1; ++i) {
    current = fabs(canonical_coefs_with_h[i] - y[i]); 
    if (current > h && current > *H) {
      *H_pos = i;
      *H = current;
    }
  }
}

static void ValleePoussin(const double *x, const double* y, unsigned N, unsigned m, double* canonical_coefs_with_h)
{
  const double eps = 1e-5;
  double H = 0.;
  unsigned H_pos = 0u;
  double h = 0.;
  // m + 1 coeffs and h at the end
  double *basis = (double *)malloc((m + 2u) * sizeof(double));
  if (!basis) {
    fprintf(stderr, "Not enough memory for Vallee-Poussin algorithm\n");
    return;
  }

  do
  {
    SelectBasis(y, N, m, basis);
    FindCanonicalCoefficients(x, basis, m, canonical_coefs_with_h);
    h = canonical_coefs_with_h[m + 1];
    FindMaximumDeviationAndItsPosition(y, canonical_coefs_with_h, m, &H, &H_pos);
    // TODO: adjust basis using algorithm described by kornev
  } while (h + eps < H);

  free(basis);
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  unsigned m = 0;
  double *x = NULL;
  double *y = NULL;
  double *canonical_coefs_with_h = NULL;

  if (argc != 2)
    return Usage(argv[0], IncorrectUsage);

  if (sscanf(argv[1], "%u", &m) != 1) {
    fprinrf(stderr, "error: parsing m from args\n");
    return Usage(argv[0], InputError);
  }

  if (fscanf(stdin, "%u", &N) != 1) {
    fprintf(stderr, "error: parsing N from input\n");
    return Usage(argv[0], InputError);
  }

  if (N < 2) {
    fprintf(stderr, "error: N(=%u) should be greater than 1\n", N);
    return Usage(argv[0], InputError);
  }

  if (m >= N - 1) {
    fprintf(stderr, "error: m(=%u) should be less than N-1(=%u)\n", m, N - 1);
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

  // m + 1 coeffs and h at the end
  canonical_coefs_with_h = (double *)malloc((m + 2u) * sizeof(double));
  if (!canonical_coefs_with_h) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    free(canonical_coefs_with_h);
    return NotEnoughMemory;
  }

  if (m == N - 2)
    FindCanonicalCoefficients(x, y, N, canonical_coefs_with_h);
  else
    ValleePoussin(x, y, N, m, canonical_coefs_with_h);

  PrintResult(x, canonical_coefs_with_h, N);

  free(x);
  free(y);
  free(canonical_coefs_with_h);
  return Success;
}
