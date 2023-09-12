#include "ExactSolution.h"

#include <math.h>

double ExactSolution(double A, double x_k) { return exp(-A * x_k); }
