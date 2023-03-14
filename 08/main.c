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
    "\tCalculate the approximate value of the 2-dimensional definite\n"
    "\tintegral using the method of quadrature formulas\n"
    , argv0);
  return error;
}

int main(int argc, const char* argv[])
{
}
