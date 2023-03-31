#include "SystemOfEquations.h"

#include <math.h>

void InitialApproximation(double *x, unsigned m) {
  for (unsigned i = 0; i < m; i++) {
    x[i] = 0.;
  }
}

static inline double RightPart(double x, double y) {
	return (y + 1) * sin(M_PI * x);
}

// FIXME: figure out by what variable derivative is taken
static inline double RightPartDerivative(double x, double y) {
	(void)y;
	return sin(M_PI * x);
}

const double X = 1.;
const double a = 0.;
const double b = 0.;

void F(double* y, const double* x, unsigned m) {
  double xi = 0.;
  const double h = X / m;

  y[0] = a;
  y[m - 1] = b;
  for (unsigned i = 1; i < m - 1; i++) {
	xi = i * h;
	y[i] = (x[i + 1] - 2. * x[i] + x[i - 1]) / h / h - RightPart(xi, x[i]);
  }
}

void dF(double* jacobian, const double* x, unsigned m) {
  double xk = 0.;
  double diff_scheme_derivative = 0.;
  const double h = X / m;

  // FIXME: derivative is wrong, add initial conditions
  for (unsigned i = 0; i < m; i++) {
	xk = i * h;
	for (unsigned j = 0; j < m; j++) {
		if (i == j)
			diff_scheme_derivative = -2.;
		else if (i == j + 1 || i + 1 == j)
			diff_scheme_derivative = 1.;
		else
			diff_scheme_derivative = 0.;
		jacobian[i * m + j] = diff_scheme_derivative / h / h - RightPartDerivative(xk, x[j]);
	}
  }
}
