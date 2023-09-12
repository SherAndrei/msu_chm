#include "SingleStep.h"

double SingleStep(int k, double h, double A, double current) {
  if (k == 0) {
    return 1.;
  }
  return current * (1. - h * A);
}
