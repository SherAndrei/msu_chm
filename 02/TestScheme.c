#include "SingleStep.h"

#include "ExactSolution.h"

#include <math.h>

double ExactSolution(double A, double x_k) {
  (void)A;
  return (55. * exp(7.* x_k) - 7. * sin(2 * x_k) - 2. * cos(2. * x_k)) / 53.;
}

double SingleStep(int k, double h, double A, double current) {
  const double x_k = k*h;
  if (k == 0) {
    return 1.;
  }
  (void)A;
  return current * (1. + 7. * h + 49. * h * h / 2.) + h / 2. * (sin(2. * x_k) + sin(2. * (x_k + h)) + 7. * h * sin(2. * x_k));
}
