#include "Error.h"
#include "GaussianElimination.h"

#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  printf(
      "Usage: %s\n"
      "\tSolve linear equations Ax=y using gaussian elimination\n"
      "\tchoosing max element by columns.\n"
      "\n"
      "\tExpected input format:\n"
      "\t\tN\n"
      "\t\ta11\ta12\t...\ta1N\ty1\n"
      "\t\ta21\ta22\t...\ta2N\ty2\n"
      "\t\t...\t...\t...\t...\t...\n"
      "\t\taN1\taN2\t...\taNN\tyN\n"
      "\n"
      "\tOutput format of solution to equation if matrix is not degenerate:\n"
      "\t\tx1\tx2\t...\txN\n",
      argv0);
  return error;
}

static int ReadSize(unsigned *N) {
  if (fscanf(stdin, "%u", N) != 1)
    return -1;
  return 0;
}

static int ReadLinearEquations(unsigned N, double *A, double *y) {
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      if (fscanf(stdin, "%lf", A + i * N + j) != 1)
        return -1;
    }
    if (fscanf(stdin, "%lf", y + i) != 1)
      return -1;
  }
  return 0;
}

static void Print(const double *arr, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    fprintf(stdout, "%e%c", arr[i], " \n"[i == N - 1]);
  }
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  double *x = NULL;
  double *y = NULL;
  double *A = NULL;

  if (argc != 1)
    return Usage(argv[0], IncorrectUsage);

  if (ReadSize(&N) < 0) {
    fprintf(stderr, "error: parsing N from input\n");
    return Usage(argv[0], InputError);
  }

  x = (double *)malloc(N * sizeof(double));
  y = (double *)malloc(N * sizeof(double));
  A = (double *)malloc(N * N * sizeof(double));
  if (!x || !y || !A) {
    fprintf(stderr, "Not enough memory\n");
    free(x);
    free(y);
    free(A);
    return NotEnoughMemory;
  }

  if (ReadLinearEquations(N, A, y) < 0) {
    fprintf(stderr, "error: parsing linear equations\n");
    free(x);
    free(y);
    free(A);
    return Usage(argv[0], InputError);
  }

  GaussMaxCol(A, N, y, x);
  Print(x, N);

  free(x);
  free(y);
  free(A);
  return Success;
}
