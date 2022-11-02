#include "step.h"

#include <cfenv>
#include <cmath>
#include <cstdio>

double En(unsigned n, double A);
double ExactSolution(double A, double x_k);
void Usage(const char* argv0);

double ExactSolution(double A, double x_k)
{
	return std::exp(-A * x_k);
}

double En(unsigned n, double A)
{
	const unsigned N = std::pow(10, n);
	const double h = std::pow(10, -static_cast<double>(n));
	const double eps = 1e-15;
	const double limit = 1e200;
	double el;
	double next;
	double max_e = 0.;

	for (auto k = 0u; k <= N; k++) {
		Step(k, h, A, el, next);

		if (std::abs(next) > limit)
			return +HUGE_VAL;

		if (std::abs(next) < eps)
			return 0.;

		double diff = std::abs(next - ExactSolution(A, (k * h)));
		if (diff > max_e)
			max_e = diff;
	}
	return max_e;
}

void Usage(const char* argv0) {
	std::printf(
		"Usage: %s n A\n"
		"\tunsigned n: 10^n of [0, 1] segment divisions\n"
		"\tdouble A: task parameter\n"
		, argv0);
}

int main(int argc, const char* argv[])
{
	if (argc != 3) {
		Usage(argv[0]);
		return 1;
	}

	unsigned n = 0;
	double A = 0.;
	if (!(std::sscanf(argv[1], "%u", &n) == 1
	   && std::sscanf(argv[2], "%lf", &A) == 1)) {
		std::printf("Error: cannot parse input parameters\n");
		return 2;
	}

	// remove FPE for std::exp for small values
	std::fenv_t env;
	std::feholdexcept(&env);

	std::printf("%8e", En(n, A));		
}
