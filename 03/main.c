#include "system_of_equations.h"
#include "step.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void Usage(const char* argv0);
FILE* OpenFile(int argc, const char* argv[]);

void Print(double x, const double* y, const double* yn, unsigned n, FILE* out);
void En(double h, FILE* out);

int main(int argc, const char* argv[])
{
	double h = 0;
	FILE* out = NULL;
	if (argc < 2 || argc > 3) {
		Usage(argv[0]);
		return 1;
	}

	if (!(sscanf(argv[1], "%lf", &h) == 1 && h > 0. && h < 1.)) {
		fprintf(stderr, "Incorrect h\n");
		return 2;
	}

	out = OpenFile(argc, argv);
	if (!out) {
		fprintf(stderr, "Cannot open output file\n");
		return 3;
	}

	En(h, out);
	fclose(out);
}

void Usage(const char* argv0) {
	fprintf(stderr,
		"Usage: %s h [filename]\n"
		"\tdouble h: initial value of step between (0, 1)\n"
		"\tfilename: output file, default -- stdout\n"
		"Result output is the table in format:\n"
		"x | yn | exact solution | error\n"
		, argv0);
}

FILE* OpenFile(int argc, const char* argv[])
{
	if (argc == 2)
		return stdout;
	return fopen(argv[2], "w");
}

void Print(double x, const double* yn, const double* y_exact, unsigned n, FILE* out)
{
	double error = 0.;
	fprintf(out, "%e\t", x);
	for (unsigned i = 0; i < n; i++) {
		error = fabs(y_exact[i] - yn[i]);
		fprintf(out, "%e\t%e\t%e%c", yn[i], y_exact[i], error, "\t\n"[i == n - 1]);
	}
}

void En(double h, FILE* out)
{
	double x = 0.;
	const double x_end = 1.;
	const unsigned n = NumberOfEquations();
	double* const f       = (double*)malloc(sizeof(double) * n);
	double* const y_exact = (double*)malloc(sizeof(double) * n);
	double* const y_prev  = (double*)malloc(sizeof(double) * n);
	double* const y_curr  = (double*)malloc(sizeof(double) * n);

	if (!y_curr || !y_prev || !f || !y_exact) {
		free(f); free(y_exact); free(y_prev); free(y_curr);
		return;
	}

	ExactSolution(y_prev, x);
	ExactSolution(y_exact, x);
	Print(x, y_prev, y_exact, n, out);

	for (x = h; x <= x_end; x += h)
	{
		Step(y_curr, y_prev, n, x, h);

		ExactSolution(y_exact, x);
		Print(x, y_curr, y_exact, n, out);

		memcpy(y_prev, y_curr, sizeof(double) * n);
	}

	free(f); free(y_exact); free(y_prev); free(y_curr);
}
