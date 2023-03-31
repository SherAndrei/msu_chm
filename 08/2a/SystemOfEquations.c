#include "SystemOfEquations.h"

#include <math.h>

void F(double* y, const double* x, unsigned m) {
  (void)m;
  y[0] = x[0] + 2.;
}

void dF(double* jacobian, const double* x, unsigned m)
{
  (void)m;
  (void)x;
  jacobian[0] = 1.;
}
