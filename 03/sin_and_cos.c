#include "system_of_equations.h"

#include <math.h>

unsigned NumberOfEquations()
{
	return 2;
}

void ExactSolution(double* y, double x)
{
	y[0] = sin(x);
	y[1] = cos(x);
}

void RightPartOfEquations(double* f, const double* y, double x)
{
	(void)x;
	f[0] =  y[1];
	f[1] = -y[0];
}
