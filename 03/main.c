#include "system_of_equations.h"
#include "step.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

	if (sscanf(argv[1], "%lf", &h) != 1) {
		perror("Cannot parse h");
		return 2;
	}

	out = OpenFile(argc, argv);
	if (!out) {
		perror("Cannot open output file");
		return 3;
	}

	En(h, out);
	fclose(out);
}

void Usage(const char* argv0) {
	fprintf(stderr,
		"Usage: %s h [filename]\n"
		"\tdouble h: initial value of step\n"
		"\tfilename: output file, default -- stdout\n"
		, argv0);
}

FILE* OpenFile(int argc, const char* argv[])
{
	if (argc == 2)
		return stdout;
	return fopen(argv[2], "w");
}

void Print(double x, const double* y, const double* yn, unsigned n, FILE* out)
{
	double error = 0.;
	fprintf(out, "%e\t", x);
	for (unsigned i = 0; i < n; i++) {
		error = fabs(y[i] - yn[i]);
		fprintf(out, "%e\t%e\t%e%c", y[i], yn[i], error, "\t\n"[i == n - 1]);
	}
}

void En(double h, FILE* out)
{
	double x = 0.;
	const double x_end = 1.;
	const unsigned n = NumberOfEquations();
	double* const y  = (double*)malloc(sizeof(double) * n);
	double* const f  = (double*)malloc(sizeof(double) * n);
	double* const yn = (double*)malloc(sizeof(double) * n);

	if (!y || !yn || !f) {
		free(y), free(yn), free(f);
		return;
	}

	while (x + h <= x_end)
	{
		ExactSolution(y, x);
		Step(yn, y, n, x, h);
		Print(x, y, yn, n, out);
		x += h;
	}

	free(y), free(yn), free(f);
}
