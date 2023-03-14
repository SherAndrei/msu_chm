#include "IntegralOnSegment.h"

double IntegralOnSegment(double a, double b, double(*f)(double))
{
  const double coeff = (b - a) / 6.;
  const double x_zero  = (a + b) / 2.;
  const double x_minus = a;
  const double x_plus  = b;
  return coeff * (f(x_minus) + 4. * f(x_zero) + f(x_plus));
}
