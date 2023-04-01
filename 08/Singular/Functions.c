#include "Functions.h"

#include <math.h>

// see Newton.pdf 1a
double f(double x) { return pow(x - 1., 2); }
double df(double x) { return 2. * (x - 1.); }

// see Newton.pdf 1b
// const double h = 1e-5;
// static double approximate_derivative(double x) { return (f(x + h) - f(x - h)) / (2. * h); }

// double f(double x) { return pow(x, 3.) * pow(cos(x), 2); }
// double df(double x) { return approximate_derivative(x); }
