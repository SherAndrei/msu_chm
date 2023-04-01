#include "SystemOfEquations.h"

#include <math.h>

#define J(i, j) jacobian[(i)*m + (j)]

void InitialApproximation(double *x, unsigned m) {
  (void)m;
  x[0] = 0;
  x[1] = 0.5;
  x[2] = 1.;
}

void F(double *f, const double *x, unsigned m) {
  (void)m;
  f[0] = pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2) - 10;
  f[1] = -pow(x[0], 2) + pow(x[1] - 3, 2) - x[2];
  f[2] = x[0] + 2. - x[1] / 3. + x[2];
}

void dF(double *jacobian, const double *x, unsigned m) {
  (void)m;
  J(0, 0) = 2. * x[0];
  J(0, 1) = 2. * x[1];
  J(0, 2) = 2. * x[2];
  J(1, 0) = -2. * x[0];
  J(1, 1) = 2. * (x[1] - 3.);
  J(1, 2) = -1.;
  J(2, 0) = 1.;
  J(2, 1) = -1. / 3.;
  J(2, 2) = 1.;
}
