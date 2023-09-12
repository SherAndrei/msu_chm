#include "DoubleStep.h"

double DoubleStep(int k, double h, double A, double prev, double current) {
  if (k == 0)
    return 1.;
  if (k == 1)
    return 1. - A * h;

  return prev - 2. * A * h * current;
}
