#include "GaussianInversion.h"

#include <math.h>
#include <stdio.h>

#define A(i, j) (A[(i)*N + (j)])
#define E(i, j) (A_inv[(i)*N + (j)])

static inline void PrintMatrixAndInversed(const double *A, const double *A_inv, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      fprintf(stderr, "%14e%c", A(i, j), "\t "[j == N - 1]);
    }
    fprintf(stderr, "| ");
    for (unsigned j = 0; j < N; j++) {
      fprintf(stderr, "%14e%c", E(i, j), "\t\n"[j == N - 1]);
    }
  }
  fprintf(stderr, "\n");
}

static inline void Swap(double *const lhs, double *const rhs) {
  double temp = *lhs;
  *lhs = *rhs;
  *rhs = temp;
}

int GaussMaxCol(double *A, double *A_inv, const unsigned N) {
  unsigned i, j, k;
  double max_elem_in_columns = 0.;
  unsigned column_with_max_elem = 0;
  const double error = 1e-15;
  double c = 0.;

  // make A_inv singular
  for (i = 0; i < N; ++i)
    for (j = 0; j < N; ++j)
      A_inv[i * N + j] = (i == j);

  for (j = 0; j < N; j++) {
    max_elem_in_columns = fabs(A(j, j));
    column_with_max_elem = j;

    // find max elem between columns
    for (i = j; i < N; i++) {
      if (fabs(A(i, j)) > max_elem_in_columns) {
        column_with_max_elem = i;
        max_elem_in_columns = fabs(A(i, j));
      }
    }

    if (max_elem_in_columns < error) {
      // matrix is degenerate!
      return -1;
    }

    // swap current row with row with max elem
    if (column_with_max_elem != j) {
      for (i = 0; i < j; ++i) {
        Swap(&(E(j, i)), &(E(column_with_max_elem, i)));
      }
      for (i = j; i < N; ++i) {
        Swap(&(A(j, i)), &(A(column_with_max_elem, i)));
        Swap(&(E(j, i)), &(E(column_with_max_elem, i)));
      }
    }

    // divide all row by first elem
    c = 1. / A(j, j);
    for (i = 0; i < j; i++) {
      E(j, i) *= c;
    }
    for (i = j; i < N; i++) {
      A(j, i) *= c;
      E(j, i) *= c;
    }

    // substract all elements under diagonal
    for (i = j + 1; i < N; i++) {
      c = A(i, j);
      for (k = 0; k < j; k++) {
        E(i, k) -= c * E(j, k);
      }
      for (k = j; k < N; k++) {
        A(i, k) -= c * A(j, k);
        E(i, k) -= c * E(j, k);
      }
    }
  }

  // reverse step
  for (j = N; j >= 1; --j) {
    for (i = 0; i < j - 1; ++i) {
      for (k = 0; k < N; ++k) {
        E(i, k) -= A(i, j - 1) * E(j - 1, k);
      }
      A(i, j - 1) = 0.;
    }
  }

  return 0;
}
