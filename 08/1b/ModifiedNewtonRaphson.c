#include "Error.h"
#include "Functions.h"

#include <math.h>
#include <stdio.h>

// FIXME: delete this file, specify h as global in 1a

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
          "Usage: %s x0 eps h [z]\n"
          "DESCRIPTION:\n"
          "\tFind a numerical solution of the equation f(x)=0 with an\n"
          "\tinitial approximation `x0` and with an accuracy of `eps` using\n"
          "\tmodified Newton's method: replacing the exact calculation of the value\n"
          "\tderivative to an approximate analog.\n"
          "\tIf root is known and is specified as parameter `z`\n"
          "\tconvergence rate diagnostics shown.\n"
          "\n"
          "OPTIONS:\n"
          "\tdouble x0 -- initial approximation of a root\n"
          "\tdouble eps -- precision\n"
          "\tdouble h -- differential width\n"
          "\tdouble z -- known root (optional parameter)\n",
          argv0);
  return error;
}

static double ApproximateDerivative(double (*f)(double), double x, double h) {
  return (f(x + h) - f(x - h)) / (2. * h);
}

static int ModifiedNewtonRaphson(double *x0, double h, double (*f)(double), double eps,
                         double *z) {
  unsigned n_steps = 0u;
  const unsigned step_limit = 200u;
  const double machine_eps = 1e-14;
  double x1;
  double f0;
  double df0;

  if (z != NULL) {
    fprintf(stdout, "step\tx0\t\tf(x0)\t\tdf(x0)\t\tfabs((x0-z)/(x1-z))\n");
  }

  while (n_steps++ < step_limit) {
    f0 = f(*x0);
    df0 = ApproximateDerivative(f, *x0, h);

    if (fabs(df0) < machine_eps) {
      fprintf(stdout, "failed: derivative is too close to 0, result'll be insufficient\n");
      return -1;
    }

    x1 = *x0 - (f0 / df0);

    if (fabs(x1 - *x0) < eps)
      return 0;

    if (z != NULL) {
      fprintf(stdout, "%u\t%e\t%e\t%e\t%e\n", n_steps, *x0, f0, df0, fabs((*x0 - *z) / (x1 - *z)));
    }

    *x0 = x1;
  }

  fprintf(stdout, "failed: step limit exceeded\n");
  return -1;
}

int main(int argc, const char *argv[]) {
  double x0;
  double eps;
  double z;
  double h;
  int error;

  if (!(argc == 4 || argc == 5))
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &x0) == 1 && sscanf(argv[2], "%lf", &eps) == 1 && sscanf(argv[3], "%lf", &h) == 1)) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return InputError;
  }

  if (argc == 5 && sscanf(argv[4], "%lf", &z) != 1) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return InputError;
  }

  error = ModifiedNewtonRaphson(&x0, h, f, eps, argc == 5 ? &z : NULL);
  if (!error)
    fprintf(stdout, "root: %14e\n", x0);
  return error;
}
