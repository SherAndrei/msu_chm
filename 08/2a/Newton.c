#include "Error.h"
#include "GaussianInversion.h"
#include "SystemOfEquations.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
          "Usage: %s eps m\n"
          "DESCRIPTION:\n"
          "\tFind a numerical vector solution of system of `m` equations \n"
          "\tF(x)=0, x=(x^1,...,x^m)^T with an accuracy of `eps` using\n"
          "\tNewton's method.\n"
          "\n"
          "OPTIONS:\n"
          "\tdouble eps > 0 -- precision\n"
          "\tunsigned m > 0 -- amount of equations in system\n",
          argv0);
  return error;
}

static void Print(const double *arr, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(stdout, "%e%c", arr[i], " \n"[i == N - 1]);
  }
}

static double UniformMetric(const double* lhs, const double* rhs, unsigned N) {
  double max = 0;
  double current = 0;
  for (unsigned i = 0; i < N; i++) {
    current = fabs(rhs[i] - lhs[i]);
    if (current > max)
      max = current;
  }
  return max;
}

static void MatrixOnVectorProduct(const double* matrix, const double* vector, unsigned N, double* result)
{
  for (unsigned i = 0; i < N; i++) {
    result[i] = 0.;
    for (unsigned j = 0; j < N; j++) {
      result[i] += matrix[i * N + j] * vector[j];
    }
  }
}

static int Root(double *x_prev, void (*F)(double *, const double *, unsigned),
                void (*dF)(double *, const double *, unsigned), unsigned m, double eps) {
  const unsigned step_limit = 200u;
  unsigned counter = 0u;
  double *x_curr = NULL;
  double *y = NULL;
  double *jacobian = NULL;
  double *inversed_jacobian = NULL;

  x_curr = malloc(m * sizeof(*x_curr));
  y = malloc(m * sizeof(*y));
  jacobian = malloc(m * m * sizeof (*jacobian));
  inversed_jacobian = malloc(m * m * sizeof (*jacobian));

  if (!x_curr || !y || !jacobian || !inversed_jacobian) {
    free(x_curr);
    free(y);
    free(jacobian);
    free(inversed_jacobian);
    return NotEnoughMemory;
  }

  while (1)
  {
    if (counter++ >= step_limit) {
	    fprintf(stdout, "Limit of steps exceeded\n");
	    break;
    }
    F(y, x_prev, m);
    dF(jacobian, x_prev, m);

    // (F'(x_n))^-1
    if (GaussMaxCol(jacobian, inversed_jacobian, m) != 0) {
      fprintf(stdout, "Jacobian matrix appears to be degenerate\n");
      break;
    }

    // (F'(x_n))^-1 * F(x_n)
    MatrixOnVectorProduct(inversed_jacobian, y, m, x_curr);

    // x_{n+1} = x_n - (F'(x_n))^-1 * F(x_n)
    for (unsigned i = 0; i < m; i++) {
      x_curr[i] = x_prev[i] - x_curr[i];
    }

    if (UniformMetric(x_prev, x_curr, m) < eps)
      break;

    for (unsigned i = 0; i < m; i++) {
      x_prev[i] = x_curr[i];
    }
  }

  free(x_curr);
  free(y);
  free(jacobian);
  free(inversed_jacobian);
  return 0;
}

int main(int argc, const char *argv[]) {
  double eps;
  unsigned m;
  double *x;
  int root_search_result;

  if (argc != 3)
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &eps) == 1 && sscanf(argv[2], "%u", &m) == 1)) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return Usage(argv[0], InputError);
  }

  if (m == 0 || eps <= 0) {
    fprintf(stderr, "error: input should be greater than 0\n");
    return Usage(argv[0], InputError);
  }

  x = malloc(m * sizeof(*x));
  if (!x) {
    fprintf(stderr, "error: not enough memory\n");
    return NotEnoughMemory;
  }

  InitialApproximation(x, m);
  root_search_result = Root(x, F, dF, m, eps);
  if (root_search_result == 0)
    Print(x, m);

  free(x);
  return root_search_result;
}
