#include "Error.h"
#include "GaussianElimination.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int Usage(const char* argv0, int error) {
  printf(
      "Usage: %s degree\n"
      "DESCRIPTION:\n"
      "\tCalculate interpolation polynom of degree `degree` < (N-1) by using input data which consists N > 1 points.\n"
      "OPTIONS:\n"
      "\tunsigned degree -- degree of desired interpolation polynom\n"
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

static void FillVandermondeMatrixAndColumnWithH(const double* x, unsigned N, double* A) {
  for (unsigned i = 0u; i < N; i++) {
    for (unsigned j = 0u; j < N - 1; j++) {
      A[i * N + j] = pow(x[i], j);
    }
    A[i * N + N - 1] = pow(-1., i);
  }
}

static int FindCanonicalCoefficients(const double* x, const double* basis, unsigned n_coeffs_with_h, double* a) {
  // last column filled up with (-1)^i
  double* A = (double*)malloc(n_coeffs_with_h * n_coeffs_with_h * sizeof(double));
  if (!A) {
    fprintf(stderr, "Not enough memory for matrix\n");
    return 4;
  }

  FillVandermondeMatrixAndColumnWithH(x, n_coeffs_with_h, A);

  if (GaussMaxCol(A, n_coeffs_with_h, basis, a) < 0) {
    fprintf(stderr, "Matrix is degenerate!\n");
    free(A);
    return 1;
  }

  free(A);
  return 0;
}

static double CanonicalForm(const double* a, double x, unsigned n_coeffs_with_h) {
  double res = 0.;
  for (unsigned i = 0; i < n_coeffs_with_h; i++)
    res += a[i] * pow(x, i);
  return res;
}

static void PrintResult(const double* x, const double* y, const double* coeffs_with_h, unsigned n_coeffs_with_h) {
  double canonical_form = 0.;
  for (unsigned i = 0; i < n_coeffs_with_h; i++) {
    canonical_form = CanonicalForm(coeffs_with_h, x[i], n_coeffs_with_h);
    fprintf(stdout, "%20e %20e %20e %20e\n", x[i], y[i], canonical_form, fabs(canonical_form - y[i]));
  }
  fprintf(stdout, "\nh = %e\n", coeffs_with_h[n_coeffs_with_h - 1]);
}

// from N points select degree+1 equally distributed points for basis
static void SelectBasis(const double* y, unsigned N, unsigned n_coeffs_with_h, double* basis, unsigned* basis_indices) {
  unsigned step = N / (n_coeffs_with_h - 1u);
  unsigned i = 0u;
  unsigned basis_i = 0u;
  for (i = 0u, basis_i = 0u; i < N && basis_i < n_coeffs_with_h; i += step, ++basis_i) {
    basis[basis_i] = y[i];
    basis_indices[basis_i] = i;
  }
  assert(basis_i == n_coeffs_with_h - 1);
}

static double Delta(unsigned i, const double* x, const double* y, const double* coeffs_with_h, unsigned n_coeffs_with_h)
{
  return y[i] - CanonicalForm(coeffs_with_h, x[i], n_coeffs_with_h);
}

static void AdjustBasis(const double* y, const double* x, const double* coeffs_with_h, unsigned n_coeffs, unsigned max_deviation_pos, double* basis,
			unsigned* basis_indices) {
  if (max_deviation_pos < basis_indices[0]) {
    if (!!signbit(Delta(basis_indices[0], x, y, coeffs_with_h, n_coeffs + 1)) == !!signbit(Delta(max_deviation_pos, x, y, coeffs_with_h, n_coeffs + 1)))
      basis_indices[0] = max_deviation_pos;
    else
      return;
      // TODO: shift left
    return;
  }
  if (max_deviation_pos > basis_indices[n_coeffs - 1]) {
    return;
  }

  (void)y;
  (void)n_coeffs;
  (void)n_coeffs;
  (void)basis;
}

// if maximum deviation is less than h, then H and H_pos is unchanged
static void FindMaximumDeviationAndItsPosition(const double* y, const double* x, const double* coeffs_with_h, unsigned n_coeffs_with_h,
					       double* H, unsigned* H_pos) {
  const double h = coeffs_with_h[n_coeffs_with_h - 1];
  double current = 0.;
  for (unsigned i = 0; i < n_coeffs_with_h - 1; ++i) {
    current = fabs(Delta(i, x, y, coeffs_with_h, n_coeffs_with_h));
    if (current > h && current > *H) {
      *H_pos = i;
      *H = current;
    }
  }
}

static void ValleePoussin(const double* x, const double* y, unsigned N, unsigned n_coeffs_with_h, double* coeffs_with_h) {
  const double eps = 1e-5;
  // current deviation between current interpolation polynom and exact solution
  double h = 0.;
  // maximal deviation between current interpolation polynom and exact solution
  double H = 0.;
  unsigned H_pos = 0u;

  unsigned* basis_indices = (unsigned*)malloc(n_coeffs_with_h * sizeof(unsigned));
  double* basis = (double*)malloc(n_coeffs_with_h * sizeof(double));
  if (!basis || !basis_indices) {
    fprintf(stderr, "Not enough memory for Vallee-Poussin algorithm\n");
    return;
  }

  SelectBasis(y, N, n_coeffs_with_h, basis, basis_indices);
  do {
    FindCanonicalCoefficients(x, basis, n_coeffs_with_h, coeffs_with_h);
    h = coeffs_with_h[n_coeffs_with_h - 1];
    FindMaximumDeviationAndItsPosition(y, x, coeffs_with_h, n_coeffs_with_h, &H, &H_pos);
	
	  if (fabs(H - h) < eps)
      break;

    AdjustBasis(y, x, coeffs_with_h, n_coeffs_with_h, H_pos, basis, basis_indices);
  } while (1);

  free(basis);
  free(basis_indices);
}

int main(int argc, const char* argv[]) {
  unsigned N = 0;
  unsigned polynom_degree = 0;
  double* x = NULL;
  double* y = NULL;
  // amount of coeffitients for polynom plus position for h
  unsigned n_coeffs_with_h = 0u;
  double* coeffs_with_h = NULL;

  if (argc != 2)
    return Usage(argv[0], IncorrectUsage);

  if (sscanf(argv[1], "%u", &polynom_degree) != 1) {
    fprintf(stderr, "error: parsing degree from args\n");
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

  if (polynom_degree >= N - 1) {
    fprintf(stderr, "error: degree(=%u) should be less than N-1(=%u)\n", polynom_degree, N - 1);
    return Usage(argv[0], InputError);
  }

  x = (double*)malloc(N * sizeof(double));
  y = (double*)malloc(N * sizeof(double));
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

  n_coeffs_with_h = polynom_degree + 2;
  // degree + 1 coeffs and h at the end
  coeffs_with_h = (double*)malloc(n_coeffs_with_h * sizeof(double));
  if (!coeffs_with_h) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    free(coeffs_with_h);
    return NotEnoughMemory;
  }

  if (n_coeffs_with_h == N)
    FindCanonicalCoefficients(x, y, N, coeffs_with_h);
  else
    ValleePoussin(x, y, N, polynom_degree, coeffs_with_h);

  PrintResult(x, y, coeffs_with_h, n_coeffs_with_h);

  free(x);
  free(y);
  free(coeffs_with_h);
  return Success;
}
