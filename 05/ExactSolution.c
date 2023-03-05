#include "ExactSolution.h"

#include <math.h>

// double ExactSolution(double x) { return 1 / (25 * (x * x) + 1); }
// double ExactSolution(double x) { return 1+2*x+3*pow(x,2)+4*pow(x,
// 3)+5*pow(x,4); }
double ExactSolution(double x) { return pow(x, 5); }
