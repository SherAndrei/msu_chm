#include "Error.h"
#include "Integral.h"

#include <math.h>
#include <stdio.h>

static double Integral(double a, double b, double(*f)(double), unsigned N)
{
  const double step = (b - a) / (N - 1.);
  double sum = 0.;
  for (unsigned i = 0; i < N; ++i) {
    sum += IntegralOnSegment(a, b, f);
    a += step;
    b += step;
  }
  return sum;
}

static inline double NthPower(double x)
{
  const unsigned N = 4;
  return pow(x, N);
}

static int Usage(const char* argv0, int error)
{
  fprintf(stdout,
    "Usage: %s a b N\n"
    "DESCRIPTION:\n"
    "\tCalculate the approximate value of the 1-dimensional definite\n"
    "\tintegral using the method of quadrature formulas\n"
    "\n"
    "OPTIONS:\n"
    "\tdouble a -- left bound of desired segment\n"
    "\tdouble b (b > a)-- right bound of desired segment\n"
    "\tunsigned N -- number of partitions of the segment into equal subsegments\n"
    , argv0);
  return error;
}

int main(int argc, const char* argv[])
{
  double a = 0.;
  double b = 0.;
  unsigned N = 0u;

  if (argc != 4)
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &a) == 1
     && sscanf(argv[2], "%lf", &b) == 1
     && sscanf(argv[3], "%u",  &N) == 1
     && b > a)) {
    fprintf(stderr, "error: parsing input parameters\n");
    return Usage(argv[0], InputError);
  }

  printf("%20e\n", Integral(a, b, sin, N));
}
