#include "Error.h"
#include "Integral.h"

#include <math.h>
#include <stdio.h>

static inline double Cos100(double x) { return cos(100. * x); }

static inline double ExpMinus1000(double x) { return exp(x * (-1000.)); }

static inline double Density(double x) { return 1. / sqrt(1. - x * x); }

static inline double Polynom(double x) { return 6. * pow(x, 5); }

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
          "Usage: %s a b N\n"
          "DESCRIPTION:\n"
          "\tCalculate the approximate value of the 1-dimensional definite\n"
          "\tintegral using the method of quadrature formulas\n"
          "\n"
          "OPTIONS:\n"
          "\tdouble a -- left bound of desired segment\n"
          "\tdouble b (b > a)-- right bound of desired segment\n"
          "\tunsigned N (N > 0) -- number of partitions of the segment into equal subsegments\n",
          argv0);
  return error;
}

int main(int argc, const char *argv[]) {
  double a = 0.;
  double b = 0.;
  unsigned N = 0u;
  double result = 0.;

  if (argc != 4)
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &a) == 1 && sscanf(argv[2], "%lf", &b) == 1 &&
        sscanf(argv[3], "%u", &N) == 1 && b > a && N > 0u)) {
    fprintf(stderr, "error: parsing input parameters\n");
    return Usage(argv[0], InputError);
  }

  result = Integral(a, b, Polynom, N);
  printf("%.14lf\n", result);
}
