#include "step.h"
#include "system_of_equations.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void Usage(const char *argv0) {
  fprintf(stderr,
          "Usage: %s h [filename]\n"
          "\tdouble h: initial value of step between (0, 1)\n"
          "\tfilename: output file, default -- stdout\n"
          "Result output is the table in format:\n"
          "x | yn | exact solution | error\n",
          argv0);
}

static FILE *OpenFile(int argc, const char *argv[]) {
  if (argc == 2)
    return stdout;
  return fopen(argv[2], "w");
}

static void Print(double x, const double *yn, const double *y_exact, unsigned n,
                  FILE *out) {
  fprintf(out, "%e\t", x);
  for (unsigned i = 0; i < n; i++)
    fprintf(out, "%e\t%e\t%e%c", yn[i], y_exact[i], fabs(y_exact[i] - yn[i]),
            "\t\n"[i == n - 1]);
}

static double ErrorUsingRungeRule(const double *a, const double *b,
                                  unsigned n) {
  const double runge_kutte_order = 3.;
  double sum = 0.;
  for (unsigned i = 0u; i < n; i++)
    sum += (a[i] - b[i]) * (a[i] - b[i]);
  return sqrt(sum) / (pow(2., runge_kutte_order) - 1);
}

static void En(double current_h, FILE *out) {
  const double eps_min = 1e+3;
  const double eps_max = 1e+5;
  double x = 0.;
  double next_h = 0.;
  const unsigned limit_of_steps = 10000u;
  unsigned count = 0u;
  const double end = 1.;
  const unsigned n = NumberOfEquations();
  double *const f = malloc(sizeof(*f) * n);
  double *const y_exact = malloc(sizeof(*y_exact) * n);
  double *const y_prev = malloc(sizeof(*y_prev) * n);
  double *const y_curr = malloc(sizeof(*y_curr) * n);
  double *const divided_y_curr_middle = malloc(sizeof(*divided_y_curr_middle) * n);
  double *const divided_y_curr_finished = malloc(sizeof(*divided_y_curr_finished) * n);

  if (!y_curr || !y_prev || !f || !y_exact || !divided_y_curr_middle ||
      !divided_y_curr_finished) {
    fprintf(stderr, "Not enough memory\n");
    free(f);
    free(y_exact);
    free(y_prev);
    free(y_curr);
    free(divided_y_curr_middle);
    free(divided_y_curr_finished);
    return;
  }

  ExactSolution(y_prev, x);
  ExactSolution(y_exact, x);
  Print(x, y_prev, y_exact, n, out);
  while (x < end) {
    if (count++ > limit_of_steps) {
      fprintf(stderr, "limit of steps exceeded\n");
      break;
    }

    if (x + current_h > end)
      current_h = end - x;

    Step(y_curr, y_prev, n, x, current_h);
    Step(divided_y_curr_middle, y_prev, n, x, current_h / 2.);
    Step(divided_y_curr_finished, divided_y_curr_middle, n, x + current_h / 2.,
         current_h / 2.);

    while (ErrorUsingRungeRule(divided_y_curr_finished, y_curr, n) > eps_max) {
      current_h /= 2.;
      Step(y_curr, y_prev, n, x, current_h);
      Step(divided_y_curr_middle, y_prev, n, x, current_h / 2.);
      Step(divided_y_curr_finished, divided_y_curr_middle, n,
           x + current_h / 2., current_h / 2.);
    }

    if (ErrorUsingRungeRule(divided_y_curr_finished, y_curr, n) < eps_min)
      next_h = current_h * 2.;

    ExactSolution(y_exact, x + current_h);
    Print(x + current_h, y_curr, y_exact, n, out);
    x += current_h;
    current_h = next_h;
    for (unsigned i = 0u; i < n; ++i)
      y_prev[i] = y_curr[i];
  }

  free(f);
  free(y_exact);
  free(y_prev);
  free(y_curr);
  free(divided_y_curr_middle);
  free(divided_y_curr_finished);
}

int main(int argc, const char *argv[]) {
  double h = 0;
  FILE *out = NULL;
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
