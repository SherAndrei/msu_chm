#include "step.h"
#include "system_of_equations.h"

#include <stdlib.h>

static void fill_k_impl(unsigned n, double h, const double *f, double *k) {
  for (unsigned i = 0u; i < n; i++)
    k[i] = h * f[i];
}

static void fill_k1(double *f, const double *yn, double x, double h, unsigned n,
                    double *k1) {
  RightPartOfEquations(f, yn, x);
  fill_k_impl(n, h, f, k1);
}

static void fill_k2(double *f, double *y_tmp, const double *yn, double x,
                    double h, unsigned n, const double *k1, double *k2) {
  for (unsigned i = 0u; i < n; i++)
    y_tmp[i] = yn[i] + k1[i] / 2.;
  RightPartOfEquations(f, y_tmp, x + h / 2.);
  fill_k_impl(n, h, f, k2);
}

static void fill_k3(double *f, double *y_tmp, const double *yn, double x,
                    double h, unsigned n, const double *k2, double *k3) {
  for (unsigned i = 0u; i < n; i++)
    y_tmp[i] = yn[i] + (3. * k2[i] / 4.);
  RightPartOfEquations(f, y_tmp, x + (3. * h / 4.));
  fill_k_impl(n, h, f, k3);
}

void Step(double *y_curr, const double *y_prev, unsigned n, double x,
          double h) {
  double *const f = (double *)malloc(n * sizeof(double));
  double *const y_tmp = (double *)malloc(n * sizeof(double));
  double *const k1 = (double *)malloc(n * sizeof(double));
  double *const k2 = (double *)malloc(n * sizeof(double));
  double *const k3 = (double *)malloc(n * sizeof(double));

  if (!f || !y_tmp || !k1 || !k2 || !k3) {
    free(f), free(y_tmp), free(k1), free(k2), free(k3);
    return;
  }

  fill_k1(f, y_prev, x, h, n, k1);
  fill_k2(f, y_tmp, y_prev, x, h, n, k1, k2);
  fill_k3(f, y_tmp, y_prev, x, h, n, k2, k3);

  for (unsigned i = 0u; i < n; i++) {
    y_curr[i] = y_prev[i] + (2. * k1[i] + 3. * k2[i] + 4. * k3[i]) / 9.;
  }

  free(f), free(y_tmp), free(k1), free(k2), free(k3);
}
