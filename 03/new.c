#include "system_of_equations.h"

// ./new.out 0.2
// eps_min = 1e+3
// eps_max = 1e+5

unsigned NumberOfEquations(void) { return 3; }

void ExactSolution(double *y, double x) {
  y[0] = 1 + x;
  y[1] = x * x;
  y[2] = x * x * x;
}

void RightPartOfEquations(double *f, const double *y, double x) {
  (void)y;
  f[0] = 1;
  f[1] = 2 * x;
  f[2] = 3 * x * x;
}
