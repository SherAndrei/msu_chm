#include "Task.h"

#include <math.h>

static double ExactSolutionValue(double x)
{
	return M_PI*x + sin(M_PI*x);
}

void ExactSolution(double h, double* y, unsigned N)
{
	for (unsigned i = 1u; i <= N-1; i++)
		y[i] = ExactSolutionValue(i * h);
}

static double AddendumValue(double x)
{
	(void)x;
	return 1.;
}

void Addendum(double h, double* p, unsigned N)
{
	for (unsigned i = 1u; i <= N-1; i++)
		p[i] = AddendumValue(i * h);
}

static double RightPartValue(double x)
{
	return M_PI*M_PI*sin(M_PI*x) + AddendumValue(x)*ExactSolutionValue(x);
}

void RightPart(double h, double* f, unsigned N)
{
	for (unsigned i = 1u; i <= N-1; i++)
		f[i] = RightPartValue(i * h);
}
