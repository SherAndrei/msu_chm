#include "DoubleStep.h"

double DoubleStep(int k, double h, double A, double prev, double current) {
  if (k == 0)
    return 1.;
  if (k == 1)
    return 1. - A * h;

  (void)current;
  return prev * (1 - A * h) / (1 + A * h);
}
