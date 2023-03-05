#include "system_of_equations.h"

#include <math.h>

unsigned NumberOfEquations(void) { return 1; }

void ExactSolution(double *y, double x) { y[0] = x + x * x + x * x * x; }

void RightPartOfEquations(double *f, const double *y, double x) {
  (void)y;
  f[0] = 1 + 2 * x + 3 * x * x;
}
