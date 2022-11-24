#include "system_of_equations.h"

#include <math.h>

unsigned NumberOfEquations()
{
	return 1;
}

void ExactSolution(double* y, double x)
{
	y[0] = -exp(-x);
}

void RightPartOfEquations(double* f, const double* y, double x)
{
	(void)x;
	f[0] = -1. * y[0];
}
