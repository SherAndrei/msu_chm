#include "Solve.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static double Scalar(const double* left, const double* right, unsigned N, double h)
{
	double result = 0.;
	for (unsigned i = 1u; i <= N-1; i++)
		result += left[i] * right[i] * h;
	return result;
}

static double EigenValue(const double* p, unsigned m, double h, unsigned N)
{
	const double angle = M_PI * (2. * m - 1.) / (2 * (2. * N - 1.));
	return pow(2. * sin(angle) / h, 2.) + p[m];
}

static void EigenVector(unsigned m, double* em, unsigned N)
{
	const double normalization_constant = sqrt(2.);
	const double fraction_in_angle = M_PI*(2.*m-1.)/(2.*N-1.);
	for (unsigned k = 1u; k <= N-1; ++k) {
		em[k] = normalization_constant * sin(fraction_in_angle * k);
	}
}

void Solve(const double* P, const double* F, double* Yn, double h, unsigned N)
{
	double dm = 0.;
	double cm = 0.;
	double* const em = (double*)calloc(N, sizeof(double));
	if (!em) {
		fprintf(stderr, "Not enough memory for algorithm\n");
		return;
	}

	for (unsigned m = 1u; m <= N-1; m++) 
	{
		EigenVector(m, em, N);
		dm = Scalar(F, em, N, h);
		cm = dm / EigenValue(P, m, h, N);
		for (unsigned i = 1u; i <= N-1; i++)
			Yn[i] += cm * em[i];
	}
	
	free(em);
}
