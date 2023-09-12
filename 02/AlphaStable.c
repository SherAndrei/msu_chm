#include "DoubleStep.h"

double DoubleStep(int k, double h, double A, double prev, double current) {
  if (k == 0)
    return 1.;
  if (k == 1)
    return 1.0 - A * h;
  
  return (4. * current - prev) / (3. + 2. * h * A);
}
