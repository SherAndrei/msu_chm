#include "step.h"

#include <cmath>
#include <cstdio>
#include <limits>

#ifndef m
#	error Please specify m
#endif

double E_n(unsigned n, double A);
double ExactSolution(double A, double x_k);
void Usage(const char* argv0);

double ExactSolution(double A, double x_k)
{
	return std::exp(-A * x_k);
}

double E_n(unsigned n, double A)
{
	const unsigned N = std::pow(10, n);
	const double h = std::pow(10, -static_cast<double>(n));
	const double limit = 1e200;
	double el;
	double next;
	double max_e = 0.;

	for (auto k = 0u; k <= N; k++) {
		Step(k, h, A, el, next);

		if (std::abs(next) > limit) {
			max_e = -1.;
			break;
		}
		if (std::abs(next) < std::numeric_limits<double>::epsilon()) {
			max_e = 0.;
			break;
		}
		double diff = std::abs(next - ExactSolution(A, (k * h)));
		if (diff > max_e)
			max_e = diff;
	}
	return max_e;
}

void Usage(const char* argv0) {
	std::printf(
		"Usage: %s A\n"
		"\tdouble A - task parameter\n"
		, argv0);
}

int main(int argc, const char* argv[])
{
	if (argc != 2) {
		Usage(argv[0]);
		return 1;
	}

	double A = 0.;
	if (std::sscanf(argv[1], "%lf", &A) != 1) {
		std::printf("Error: cannot parse input parameters\n");
		return 2;
	}

	std::printf("E_1 = %e\t", E_n(1, A));
	std::printf("E_2 = %e\t", E_n(2, A));
	std::printf("E_3 = %e\t", E_n(3, A));
	std::printf("E_6 = %e\t", E_n(6, A));
	std::printf("m = %u\t", m);
	std::printf("A = %lf\n", A);
}
