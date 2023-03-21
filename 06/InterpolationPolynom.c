#include "Error.h"
#include "GaussianElimination.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int Usage(const char *argv0, int error) {
  printf("Usage: %s degree\n"
         "DESCRIPTION:\n"
         "\tCalculate interpolation polynom of degree `degree` < (N-1) by "
         "using input data which consists N > 1 points.\n"
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
         "\tN rows, each constists starting nodes, exact solution, result of "
         "the interpolation polynom and their difference\n"
         "\t\tx1 \ty1 \tcanonical1\tdelta1\n"
         "\t\tx2 \ty2 \tcanonical2\tdelta2\n"
         "\t\t...\t...\t...\t...\n"
         "\t\txN \tyN \tcanonicalN \tdeltaN\n",
         argv0);
  return error;
}

// last column filled up with (-1)^i
static void FillVandermondeMatrixAndColumnWithH(const double *x, unsigned N, double *A) {
  for (unsigned i = 0u; i < N; i++) {
    for (unsigned j = 0u; j < N - 1; j++) {
      A[i * N + j] = pow(x[i], j);
    }
    A[i * N + N - 1] = pow(-1., i);
  }
}

static int FindCanonicalCoefficients(const double *x, const double *y, unsigned n_coeffs_with_h,
                                     double *a) {
  double *A = malloc(n_coeffs_with_h * n_coeffs_with_h * sizeof(*A));
  if (!A) {
    fprintf(stderr, "Not enough memory for matrix\n");
    return 4;
  }

  FillVandermondeMatrixAndColumnWithH(x, n_coeffs_with_h, A);

  if (GaussMaxCol(A, n_coeffs_with_h, y, a) < 0) {
    fprintf(stderr, "Matrix is degenerate!\n");
    free(A);
    return 1;
  }

  free(A);
  return 0;
}

static double CanonicalForm(const double *a, double x, unsigned n_coeffs) {
  // Horner's method
  double res = a[n_coeffs - 1] * x;
  for (unsigned i = n_coeffs - 2; i > 0; --i) {
    res += a[i];
    res *= x;
  }
  return res + a[0];
}

static void PrintResult(const double *x, const double *y, unsigned N, const double *coeffs_with_h,
                        unsigned n_coeffs_with_h) {
  double canonical_form = 0.;
  for (unsigned i = 0; i < N; i++) {
    canonical_form = CanonicalForm(coeffs_with_h, x[i], n_coeffs_with_h - 1);
    fprintf(stdout, "%20e %20e %20e %20e\n", x[i], y[i], canonical_form,
            fabs(canonical_form - y[i]));
  }
  fprintf(stdout, "\nh = %e\n", coeffs_with_h[n_coeffs_with_h - 1]);
}

// from N points select polynom_degree + 2 points for basis
static void SelectIndicesFromXForBasis(unsigned N, unsigned n_coeffs_with_h,
                                       unsigned *basis_indices_from_x) {
  unsigned i = 0u;
  (void)N;
  for (i = 0u; i < n_coeffs_with_h; i++) {
    basis_indices_from_x[i] = i;
  }
}

static void PrepareBasisUsingIndices(const double *x, const double *y,
                                     const unsigned *basis_indices_from_x, unsigned n_coeffs_with_h,
                                     double *basis, double *y_for_basis) {
  unsigned current_index;
  for (unsigned i = 0u; i < n_coeffs_with_h; ++i) {
    current_index = basis_indices_from_x[i];
    basis[i] = x[current_index];
    y_for_basis[i] = y[current_index];
  }
}

static double Delta(unsigned i, const double *x, const double *y, const double *coeffs_with_h,
                    unsigned n_coeffs_with_h) {
  return y[i] - CanonicalForm(coeffs_with_h, x[i], n_coeffs_with_h - 1);
}

// if direction < 0 shift left by one
// if direction > 0 shift right by one
// else do nothing
static void ShiftByOne(unsigned *array, unsigned size, int direction) {
  memmove(array + (direction > 0), array + (direction < 0), size * sizeof(*array) - sizeof(*array));
}

// return a pointer to the first element in pre-sorted arr which compares
// not less than val
static unsigned *FindFirstNotLess(unsigned *arr, unsigned size, unsigned val) {
  unsigned i = 0u;
  for (; i < size; ++i) {
    if (arr[i] >= val)
      break;
  }
  return arr + i;
}

// 1 if negative, 0 otherwise
static inline int Sign(double value) { return !!signbit(value); }

static int TryToAdjustBasisIndices(const double *y, const double *x, const double *coeffs_with_h,
                                   unsigned n_coeffs_with_h, unsigned max_deviation_pos,
                                   unsigned *basis_indices) {
  const unsigned basis_left_bound_i = basis_indices[0];
  const unsigned basis_right_bound_i = basis_indices[n_coeffs_with_h - 1];
  const int sign_in_max_deviation =
      Sign(Delta(max_deviation_pos, x, y, coeffs_with_h, n_coeffs_with_h));
  int sign;
  // used only if basis_left_bound < x_at_max_deviation < basis_right_bound
  unsigned *pointer_to_right_neighbour;
  unsigned *pointer_to_left_neighbour;

  if (max_deviation_pos < basis_left_bound_i) {
    sign = Sign(Delta(basis_left_bound_i, x, y, coeffs_with_h, n_coeffs_with_h));
    if (sign_in_max_deviation != sign)
      ShiftByOne(basis_indices, n_coeffs_with_h, 1);
    basis_indices[0] = max_deviation_pos;
    return 0;
  }

  if (max_deviation_pos > basis_right_bound_i) {
    sign = Sign(Delta(basis_right_bound_i, x, y, coeffs_with_h, n_coeffs_with_h));
    if (sign_in_max_deviation != sign)
      ShiftByOne(basis_indices, n_coeffs_with_h, -1);
    basis_indices[n_coeffs_with_h - 1] = max_deviation_pos;
    return 0;
  }

  pointer_to_right_neighbour = FindFirstNotLess(basis_indices, n_coeffs_with_h, max_deviation_pos);

  // is max_deviation_pos already in basis
  if (*pointer_to_right_neighbour == max_deviation_pos) {
    return 1;
  }

  pointer_to_left_neighbour = pointer_to_right_neighbour - 1;
  sign = Sign(Delta(*pointer_to_left_neighbour, x, y, coeffs_with_h, n_coeffs_with_h));

  if (sign_in_max_deviation == sign)
    *pointer_to_left_neighbour = max_deviation_pos;
  else
    *pointer_to_right_neighbour = max_deviation_pos;
  return 0;
}

// if maximum deviation is less than h, then H and H_pos is unchanged
static void FindMaximumDeviationAndItsPosition(const double *y, const double *x, unsigned N,
                                               const double *coeffs_with_h,
                                               unsigned n_coeffs_with_h, double *H,
                                               unsigned *H_pos) {
  double current;
  *H_pos = 0u;
  *H = fabs(Delta(0, x, y, coeffs_with_h, n_coeffs_with_h));
  for (unsigned i = 1; i < N; ++i) {
    current = fabs(Delta(i, x, y, coeffs_with_h, n_coeffs_with_h));
    if (current > *H) {
      *H_pos = i;
      *H = current;
    }
  }
}

static void ValleePoussin(const double *x, const double *y, unsigned N, unsigned polynom_degree,
                          double *coeffs_with_h) {
  const unsigned n_coeffs_with_h = polynom_degree + 2;
  const double eps = 1e-5;
  // current deviation between current interpolation polynom and exact solution
  double h = 0.;
  // maximal deviation between current interpolation polynom and exact solution
  double H = 0.;
  unsigned H_pos = 0u;

  double *basis = malloc(n_coeffs_with_h * sizeof(*basis));
  double *y_for_basis = malloc(n_coeffs_with_h * sizeof(*y_for_basis));
  unsigned *basis_indices_from_x = malloc(n_coeffs_with_h * sizeof(*basis_indices_from_x));
  if (!basis || !y_for_basis || !basis_indices_from_x) {
    free(basis);
    free(y_for_basis);
    free(basis_indices_from_x);
    fprintf(stderr, "Not enough memory for Vallee-Poussin algorithm\n");
    return;
  }

  SelectIndicesFromXForBasis(N, n_coeffs_with_h, basis_indices_from_x);
  do {
    PrepareBasisUsingIndices(x, y, basis_indices_from_x, n_coeffs_with_h, basis, y_for_basis);

    if (FindCanonicalCoefficients(basis, y_for_basis, n_coeffs_with_h, coeffs_with_h) != 0)
      break;

    FindMaximumDeviationAndItsPosition(y, x, N, coeffs_with_h, n_coeffs_with_h, &H, &H_pos);

    h = fabs(coeffs_with_h[n_coeffs_with_h - 1]);
    if (fabs(H - h) < eps)
      break;

    if (TryToAdjustBasisIndices(y, x, coeffs_with_h, n_coeffs_with_h, H_pos,
                                basis_indices_from_x) != 0)
      break;

  } while (1);

  free(basis);
  free(y_for_basis);
  free(basis_indices_from_x);
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  unsigned polynom_degree = 0;
  double *x = NULL;
  double *y = NULL;
  // amount of coeffitients for polynom plus position for h
  unsigned n_coeffs_with_h = 0u;
  double *coeffs_with_h = NULL;

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

  n_coeffs_with_h = polynom_degree + 2;
  // degree + 1 coeffs and h at the end
  coeffs_with_h = malloc(n_coeffs_with_h * sizeof(*coeffs_with_h));
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

  PrintResult(x, y, N, coeffs_with_h, n_coeffs_with_h);

  free(x);
  free(y);
  free(coeffs_with_h);
  return Success;
}
