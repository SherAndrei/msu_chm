#include "GaussianElimination.h"

#include <math.h>
#include <stdio.h>

#define A(i, j) (A[(i)*N + (j)])

static inline void PrintLinearEquations(const double* A, const double* y, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      fprintf(stderr, "%15e\t", A[i * N + j]);
    }
    fprintf(stderr, "|\t%15e\n", y[i]);
  }
  fprintf(stderr, "\n");
}

static inline void Swap(double* const lhs, double* const rhs) {
  double temp = *lhs;
  *lhs = *rhs;
  *rhs = temp;
}

int GaussMaxCol(double* A, unsigned N, const double* y, double* x) {
  unsigned i, j, k;
  double max_elem_in_columns = 0.;
  unsigned column_with_max_elem = 0;
  const double error = 1e-15;
  double c = 0.;

  for (i = 0; i < N; i++) {
    x[i] = y[i];
  }
  PrintLinearEquations(A, x, N);
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
      for (i = j; i < N; i++) {
	Swap(&(A(j, i)), &(A(column_with_max_elem, i)));
      }
      Swap(x + j, x + column_with_max_elem);
    }

    // divide all row by first elem
    c = 1. / A(j, j);
    for (i = j; i < N; i++) {
      A(j, i) *= c;
    }
    x[j] *= c;

    // substract all elements under diagonal
    for (i = j + 1; i < N; i++) {
      c = A(i, j);
      for (k = j; k < N; k++) {
	A(i, k) -= c * A(j, k);
      }
      x[i] -= c * x[j];
    }
  }

  // reverse step
  for (j = N; j >= 1; j--) {
    for (i = 0; i < j - 1; i++) {
      x[i] -= A(i, j - 1) * x[j - 1];
      A(i, j - 1) = 0.;
    }
  }
  PrintLinearEquations(A, x, N);
  return 0;
}
