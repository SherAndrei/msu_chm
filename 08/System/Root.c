#include "Root.h"
#include "Error.h"
#include "GaussianInversion.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline void Print(FILE *out, const double *arr, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(out, "%20.14lf\n", arr[i]);
  }
  fputc('\n', out);
}

static inline void PrintMatrix(FILE *out, const double *A, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      fprintf(out, "%20.14lf%c", A[i * N + j], "\t\n"[j == N - 1]);
    }
  }
  fputc('\n', out);
}

static double UniformMetric(const double *lhs, const double *rhs, unsigned N) {
  double max = 0;
  double current = 0;
  for (unsigned i = 0; i < N; i++) {
    current = fabs(rhs[i] - lhs[i]);
    if (current > max)
      max = current;
  }
  return max;
}

static void MatrixOnVectorProduct(const double *matrix, const double *vector, unsigned N,
                                  double *result) {
  for (unsigned i = 0; i < N; i++) {
    result[i] = 0.;
    for (unsigned j = 0; j < N; j++) {
      result[i] += matrix[i * N + j] * vector[j];
    }
  }
}

int Root(double *x_prev, void (*F)(double *, const double *, unsigned),
                void (*dF)(double *, const double *, unsigned), unsigned m, double eps,
                int verbose) {
  const unsigned step_limit = 200u;
  unsigned counter = 0u;
  double *x_curr = NULL;
  double *f = NULL;
  double *jacobian = NULL;
  double *inversed_jacobian = NULL;
  int error = Success;

  x_curr = malloc(m * sizeof(*x_curr));
  f = malloc(m * sizeof(*f));
  jacobian = malloc(m * m * sizeof(*jacobian));
  inversed_jacobian = malloc(m * m * sizeof(*jacobian));

  if (!x_curr || !f || !jacobian || !inversed_jacobian) {
    free(x_curr);
    free(f);
    free(jacobian);
    free(inversed_jacobian);
    return NotEnoughMemory;
  }

  while (1) {

    if (verbose) {
      fprintf(stderr, "step: %u\n", counter);
    }

    if (counter++ >= step_limit) {
      fprintf(stdout, "Limit of steps exceeded\n");
      error = -1;
      break;
    }

    F(f, x_prev, m);

    if (verbose) {
      fprintf(stderr, "F returned:\n");
      Print(stderr, f, m);
    }

    dF(jacobian, x_prev, m);

    if (verbose) {
      fprintf(stderr, "Jacobian is:\n");
      PrintMatrix(stderr, jacobian, m);
    }

    // (F'(x_n))^-1
    if (GaussMaxCol(jacobian, inversed_jacobian, m) != 0) {
      fprintf(stdout, "Jacobian matrix appears to be degenerate\n");
      error = -1;
      break;
    }

    if (verbose) {
      fprintf(stderr, "Inversed Jacobian is:\n");
      PrintMatrix(stderr, inversed_jacobian, m);
    }

    // (F'(x_n))^-1 * F(x_n)
    MatrixOnVectorProduct(inversed_jacobian, f, m, x_curr);

    // x_{n+1} = x_n - (F'(x_n))^-1 * F(x_n)
    for (unsigned i = 0; i < m; i++) {
      x_curr[i] = x_prev[i] - x_curr[i];
    }

    if (verbose) {
      fprintf(stderr, "Solution on current step is:\n");
      Print(stderr, x_curr, m);
    }

    if (UniformMetric(x_prev, x_curr, m) < eps)
      break;

    for (unsigned i = 0; i < m; i++) {
      x_prev[i] = x_curr[i];
    }
  }

  free(x_curr);
  free(f);
  free(jacobian);
  free(inversed_jacobian);
  return error;
}
