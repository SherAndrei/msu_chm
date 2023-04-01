#include "SystemOfEquations.h"

#include <math.h>

// See Newton.pdf 2a

void InitialApproximation(double *x, unsigned m) {
  for (unsigned i = 0; i < m; i++) {
    x[i] = 0.;
  }
}

static double RightPart(double x, double y) { return cos(x) * sin(y); }

static double RightPartDerivativeByY(double x, double y) { return cos(x) * cos(y); }

const double X = 10;
const double a = -1;
const double b = 1;

// amount of algebraic equations is N - 1
// amounf of equations required for Newton is m
// m == N - 1

void F(double *y, const double *x, unsigned m) {
  const double h = X / (m + 1);
  const double divide_by_hh = 1. / h / h;

  y[0] = (x[1] - 2. * x[0] + a) * divide_by_hh - RightPart(h, x[0]);
  y[m - 1] = (b - 2. * x[m - 1] + x[m - 2]) * divide_by_hh - RightPart(X - h, x[m - 1]);
  for (unsigned k = 1; k < m - 1; k++) {
    y[k] = (x[k + 1] - 2. * x[k] + x[k - 1]) * divide_by_hh - RightPart((k + 1) * h, x[k]);
  }
}

#define J(i, j) jacobian[(i)*m + (j)]

void dF(double *jacobian, const double *x, unsigned m) {
  const double h = X / (m + 1);
  const double main_diagonal = -2. / h / h;
  const double subdiagonal = 1. / h / h;
  const double supradiagonal = subdiagonal;

  for (unsigned i = 0; i < m * m; i++) {
    jacobian[i] = 0.;
  }

  J(0, 0) = main_diagonal - RightPartDerivativeByY(h, x[0]);
  J(0, 1) = supradiagonal;
  J(m - 1, m - 1) = main_diagonal - RightPartDerivativeByY(X - h, x[m - 1]);
  J(m - 1, m - 2) = subdiagonal;
  for (unsigned k = 1; k < m - 1; k++) {
    J(k, k - 1) = subdiagonal;
    J(k, k + 0) = main_diagonal - RightPartDerivativeByY((k + 1) * h, x[k]);
    J(k, k + 1) = supradiagonal;
  }
}
