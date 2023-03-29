#include "Error.h"
#include "GaussianInversion.h"

#include <stdio.h>
#include <stdlib.h>

static int Usage(const char *argv0, int error) {
  printf("Usage: %s\n"
         "\tInvert matrix A using gaussian elimination\n"
         "\tchoosing max element by columns.\n"
         "\n"
         "\tExpected input format:\n"
         "\t\tN\n"
         "\t\ta11\ta12\t...\ta1N\n"
         "\t\ta21\ta22\t...\ta2N\n"
         "\t\t...\t...\t...\t...\n"
         "\t\taN1\taN2\t...\taNN\n"
         "\n"
         "\tOutput format is inversed matrix if it is not degenerate.\n",
         argv0);
  return error;
}

static int ReadSize(unsigned *N) {
  if (fscanf(stdin, "%u", N) != 1)
    return -1;
  return 0;
}

static int ReadMatrix(unsigned N, double *A) {
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      if (fscanf(stdin, "%lf", A + i * N + j) != 1)
        return -1;
    }
  }
  return 0;
}

static void PrintMatrix(const double *matrix, unsigned N) {
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      fprintf(stdout, "%14e%c", matrix[i * N + j], "\t\n"[j == N - 1]);
    }
  }
  fprintf(stdout, "\n");
}

int main(int argc, const char *argv[]) {
  unsigned N = 0;
  double *A = NULL;
  double *A_inv = NULL;

  if (argc != 1)
    return Usage(argv[0], IncorrectUsage);

  if (ReadSize(&N) < 0) {
    fprintf(stderr, "error: parsing N from input\n");
    return Usage(argv[0], InputError);
  }

  A = malloc(N * N * sizeof(*A));
  A_inv = malloc(N * N * sizeof(*A_inv));
  if (!A_inv || !A) {
    fprintf(stderr, "Not enough memory\n");
    free(A);
    free(A_inv);
    return NotEnoughMemory;
  }

  if (ReadMatrix(N, A) < 0) {
    fprintf(stderr, "error: parsing matrix\n");
    free(A);
    free(A_inv);
    return Usage(argv[0], InputError);
  }

  GaussMaxCol(A, A_inv, N);
  PrintMatrix(A_inv, N);

  free(A);
  free(A_inv);
}
