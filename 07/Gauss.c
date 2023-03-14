#include "IntegralOnSegment.h"

#include <math.h>

double IntegralOnSegment(double a, double b, double(*f)(double))
{
  const double coeff = (b - a) / 18.;
  const double x_zero  = (a + b) / 2.;
  const double x_minus = x_zero - (b - a) / 2. * sqrt(3./5.);
  const double x_plus  = x_zero + (b - a) / 2. * sqrt(3./5.);
  return coeff * (5. * f(x_minus) + 8. * f(x_zero) + 5. * f(x_plus));
}
