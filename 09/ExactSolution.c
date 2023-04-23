#include "ExactSolution.h"

#include <math.h>

double ExactSolution(double x1, double x2) {
  return 10. * sin(10. * (pow(x1 - 0.5, 2.) + pow(x2 - 0.5, 2.)));
}
