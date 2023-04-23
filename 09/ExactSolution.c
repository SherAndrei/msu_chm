#include "ExactSolution.h"

#include <math.h>

double ExactSolution(double x, double y) {
  return 2 * sin(10 * sqrt((pow(x - 0.5, 2.) + pow(y - 0.5, 2.))));
}
