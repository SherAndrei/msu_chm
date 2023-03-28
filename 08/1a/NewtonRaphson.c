#include "Error.h"
#include "Functions.h"

#include <math.h>
#include <stdio.h>

static int Usage(const char *argv0, int error) {
  fprintf(stdout,
          "Usage: %s x0 eps [z]\n"
          "DESCRIPTION:\n"
          "\tFind a numerical solution of the equation f(x)=0 with an\n"
          "\tinitial approximation `x0` and with an accuracy of `eps` using Newton's method.\n"
          "\tIf root is known and is specified as parameter `z`\n"
          "\tconvergence rate diagnostics shown.\n"
          "\n"
          "OPTIONS:\n"
          "\tdouble x0 -- initial approximation of a root\n"
          "\tdouble eps -- precision\n"
          "\tdouble z -- known root (optional parameter)\n",
          argv0);
  return error;
}

static int NewtonRaphson(double *x0, double (*f)(double), double (*df)(double), double eps,
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
    df0 = df(*x0);

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
  int error;

  if (!(argc == 3 || argc == 4))
    return Usage(argv[0], IncorrectUsage);

  if (!(sscanf(argv[1], "%lf", &x0) == 1 && sscanf(argv[2], "%lf", &eps) == 1)) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return InputError;
  }

  if (argc == 4 && sscanf(argv[3], "%lf", &z) != 1) {
    fprintf(stderr, "error: parsing input parameters failed\n");
    return InputError;
  }

  error = NewtonRaphson(&x0, f, df, eps, argc == 4 ? &z : NULL);
  if (!error)
    fprintf(stdout, "root: %14e\n", x0);
  return error;
}
