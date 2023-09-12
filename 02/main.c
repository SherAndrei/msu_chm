#include "Scheme.h"

#include <fenv.h>
#include <stdio.h>

static void Usage(const char *argv0) {
  fprintf(stderr, "Usage: %s N X A\n"
         "\tunsigned N: amount of segment divisions\n"
         "\tdouble X: length of the segment\n"
         "\tdouble A: task parameter\n",
         argv0);
}

int main(int argc, const char *argv[]) {
  unsigned N;
  double X;
  double A;
  fenv_t env;
  double max_error;
  if (argc != 4) {
    Usage(argv[0]);
    return 1;
  }

  if (!(sscanf(argv[1], "%u", &N)  == 1
     && sscanf(argv[2], "%lf", &X) == 1
     && sscanf(argv[3], "%lf", &A) == 1)) {
    printf("Error: cannot parse input parameters\n");
    return 2;
  }

  // remove FPE for exp for small values
  feholdexcept(&env);

  max_error = Scheme(N, X, A);
  printf("# error: %e\n", max_error);
}
