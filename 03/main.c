#include "system_of_equations.h"
#include "step.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static void Usage(const char* argv0) {
	fprintf(stderr,
		"Usage: %s h [filename]\n"
		"\tdouble h: initial value of step between (0, 1)\n"
		"\tfilename: output file, default -- stdout\n"
		"Result output is the table in format:\n"
		"x | yn | exact solution | error\n"
		, argv0);
}

static FILE* OpenFile(int argc, const char* argv[])
{
	if (argc == 2)
		return stdout;
	return fopen(argv[2], "w");
}

static void Print(double x, const double* yn, const double* y_exact, unsigned n, FILE* out)
{
	fprintf(out, "%e\t", x);
	for (unsigned i = 0; i < n; i++)
		fprintf(out, "%e\t%e\t%e%c", yn[i], y_exact[i], fabs(y_exact[i] - yn[i]), "\t\n"[i == n - 1]);
}

static double ErrorUsingRungeRule(const double* a, const double* b, unsigned n)
{
	const double runge_kutte_order = 3.;
	double sum = 0.;
	for (unsigned i = 0u; i < n; i++)
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	return sum / (pow(2., runge_kutte_order) - 1);
}

static void En(double h_started, FILE* out)
{
	const double eps_min = 1e+3;
	const double eps_max = 1e+5;
	double x = 0.;
	double h_finished = 0.;
	const unsigned limit_of_steps = 10000u;
	unsigned count = 0u;
	const double end = 1.;
	const unsigned n = NumberOfEquations();
	double* const f       = (double*)malloc(sizeof(double) * n);
	double* const y_exact = (double*)malloc(sizeof(double) * n);
	double* const y_prev  = (double*)malloc(sizeof(double) * n);
	double* const y_curr  = (double*)malloc(sizeof(double) * n);
	double* const divided_y_curr_middle = (double*)malloc(sizeof(double) * n);
	double* const divided_y_curr_finished = (double*)malloc(sizeof(double) * n);

	if (!y_curr || !y_prev || !f || !y_exact || !divided_y_curr_middle || !divided_y_curr_finished) {
		fprintf(stderr, "Not enough memory\n");
		free(f); free(y_exact); free(y_prev); free(y_curr); free(divided_y_curr_middle); free(divided_y_curr_finished);
		return;
	}

	ExactSolution(y_prev, x);
	ExactSolution(y_exact, x);
	Print(x, y_prev, y_exact, n, out);
	while (x < end)
	{
		if (count++ > limit_of_steps)
		{
			fprintf(stderr, "limit of steps exceeded\n");
			break;
		}

		if (x + h_started > end)
			h_started = end - x;

		Step(y_curr, y_prev, n, x, h_started);
		Step(divided_y_curr_middle, y_prev, n, x, h_started / 2.);
		Step(divided_y_curr_finished, divided_y_curr_middle, n, x + h_started / 2., h_started / 2.);

		while (ErrorUsingRungeRule(divided_y_curr_finished, y_curr, n) > eps_max)
		{
			h_started /= 2.;
			Step(y_curr, y_prev, n, x, h_started);
			Step(divided_y_curr_middle, y_prev, n, x, h_started / 2.);
			Step(divided_y_curr_finished, divided_y_curr_middle, n, x + h_started / 2., h_started / 2.);
		}

		if (ErrorUsingRungeRule(divided_y_curr_finished, y_curr, n) < eps_min)
			h_finished = h_started * 2.;

		ExactSolution(y_exact, x + h_started);
		Print(x + h_started, y_curr, y_exact, n, out);
		x += h_started;
		h_started = h_finished;
		memcpy(y_prev, y_curr, sizeof(double) * n);
	}

	free(f); free(y_exact); free(y_prev); free(y_curr); free(divided_y_curr_middle); free(divided_y_curr_finished);
}

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
