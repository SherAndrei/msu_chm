#include "step.h"

#include <fenv.h>
#include <math.h>
#include <stdio.h>

double En(unsigned n, double A);
double ExactSolution(double A, double x_k);
void Usage(const char* argv0);

double ExactSolution(double A, double x_k)
{
	return exp(-A * x_k);
}

double En(unsigned n, double A)
{
	const unsigned N = pow(10, n);
	const double h = pow(10, -(double)n);
	const double eps = 1e-15;
	const double limit = 1e200;
	double el;
	double next;
	double diff;
	double max_e = 0.;

	for (unsigned k = 0u; k <= N; k++) {
		Step(k, h, A, &el, &next);

		if (fabs(next) > limit)
			return +HUGE_VAL;

		if (fabs(next) < eps)
			return 0.;

		diff = fabs(next - ExactSolution(A, (k * h)));
		if (diff > max_e)
			max_e = diff;
	}
	return max_e;
}

void Usage(const char* argv0) {
	printf(
		"Usage: %s n A\n"
		"\tunsigned n: 10^n of [0, 1] segment divisions\n"
		"\tdouble A: task parameter\n"
		, argv0);
}

int main(int argc, const char* argv[])
{
	unsigned n = 0;
	double A = 0.;
	fenv_t env;
	if (argc != 3) {
		Usage(argv[0]);
		return 1;
	}

	if (!(sscanf(argv[1], "%u", &n) == 1
	   && sscanf(argv[2], "%lf", &A) == 1)) {
		printf("Error: cannot parse input parameters\n");
		return 2;
	}

	// remove FPE for exp for small values
	feholdexcept(&env);

	printf("%8e", En(n, A));
}
