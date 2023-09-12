#include "Scheme.h"

#include "DoubleStep.h"
#include "ExactSolution.h"

#include <math.h>
#include <stdio.h>

double Scheme(unsigned N, double X, double A) {
  const double h = (X - 0.) / N;
  const double eps = 1e-15;
  const double limit = 1e200;
  double prevprev = 0.;
  double prev = 0.;
  double current = 0.;
  double max_e = 0.;
  double solution = 0.;

  for (unsigned k = 0u; k <= N; k++) {
    current = DoubleStep(k, h, A, prevprev, prev);

    if (fabs(current) > limit)
      return +HUGE_VAL;

    if (fabs(current) < eps)
      return 0.;

    solution = ExactSolution(A, (k * h));
    max_e = fmax(max_e, fabs(current - solution));    
    printf("%20.15lf %20.15lf %20.15lf\n", (k * h), current, solution);
  
    prevprev = prev;
    prev = current;
  }
  return max_e;
}
