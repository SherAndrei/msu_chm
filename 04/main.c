#include "Solve.h"
#include "Task.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static FILE* OpenFile(int argc, const char* argv[])
{
	if (argc == 2)
		return stdout;
	return fopen(argv[2], "w");
}

static unsigned ParseToUnsigned(const char* param_name, const char* from, unsigned def) {
	unsigned ret = def;
	if (sscanf(from, "%u", &ret) != 1) {
		fprintf(stderr, "Error parsing %s has occured, setting default value %u\n", param_name, def);
	}
	return ret;
}

static void Usage(const char* prog_name) {
	printf(
		"Usage: %s N [output.txt]\n"
		"\tunsigned N - matrix dim, N > 2\n"
		"\tfilename: output file, default -- stdout\n"
		"Solving difference equation -y''(x)+P(x)y(x)=F(x)\n"
		"with initial conditions y(0)=0, y'(1)=0 using next\n"
		"difference scheme -\\frac{y_{k-1}-2y_k+y_{k+1}}{h^2}+P_ky_k=F_k\n"
		"and corresponding initial conditions y_0=0, y_N=y_{N-1}.\n"
		"Result output is the table in format:\n"
		"x | yn | y \n"
		, prog_name);
}

// ||u_h - (u)_h||_h
static double Error(const double* solved, const double* exact, unsigned N, double h)
{
	double result = 0.;
	for (unsigned i = 0u; i < N; i++)
		result += pow(solved[i] - exact[i], 2.) * h;
	return sqrt(result);
}

int main(int argc, const char* argv[])
{
	unsigned N = 0;
	FILE* out = NULL;
	double h = 0.;
	double* right = NULL;
	double* add = NULL;
	double* solved = NULL;
	double* exact = NULL;

	if (argc < 2 || argc > 3) {
		Usage(argv[0]);
		return 1;
	}
	
	N = ParseToUnsigned("N", argv[1], 0);
	if (N < 3)
		return 2;
	
	out = OpenFile(argc, argv);
	if (!out) {
		fprintf(stderr, "Cannot open output file\n");
		return 3;
	}

	right = (double*)calloc(N, sizeof(double));
	add = (double*)calloc(N, sizeof(double));
	solved = (double*)calloc(N, sizeof(double));
	exact = (double*)calloc(N, sizeof(double));

	if (!right || !add || !solved || !exact) {
		fprintf(stderr, "Not enough memory\n");
		fclose(out);
		free(right);
		free(add);
		free(solved);
		free(exact);
		return 4;
	}

	h = 2./(2.* N - 1.);
	
	RightPart(h, right, N);
	Addendum(h, add, N);
	Solve(add, right, solved, h, N);
	ExactSolution(h, exact, N);

	fprintf(out, "%e\t%e\t%e\n", 0., 0., 0.); // initial condition
	for (unsigned i = 1u; i <= N-1; i++) {
		fprintf(out, "%e\t%e\t%e\n", i*h, solved[i], exact[i]);
	}

	printf("Error: %e\n", Error(solved, exact, N, h));

	fclose(out);
	free(right);
	free(add);
	free(solved);
	free(exact);
}
