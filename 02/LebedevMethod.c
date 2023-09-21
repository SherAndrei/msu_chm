#include "Scheme.h"

#include "ExactSolution.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// solution for y'=-A(y-sinx)+cosx
double ExactSolution(double A, double x_k) {
  return exp(-A * x_k) + sin(x_k);
}

static double SingleStep(int k, double h, double A, double x_k, double y_k) {
  if (k == 0) {
    return 1.;
  }
  
  // explicit Euler's scheme
  return y_k  * (1 - A * h) + h * (A * sin(x_k) + cos(x_k));
}

static inline double ChebyshevNode(unsigned k, unsigned N, double left, double right) {
  const double add = 0.5 * (left + right);
  const double coef = 0.5 * (right - left);
  return add + coef * cos(M_PI * (2. * (N + 1 - k) - 1.) / (2. * N));
}

static unsigned* GenerateShuffledIndices(unsigned N) {
  unsigned* indices = malloc(N * sizeof(*indices));
  unsigned* tmp = malloc(N * sizeof(*indices));

  if (!indices || !tmp) {
    free(tmp);
    free(indices);
    return NULL;
  }

  *tmp = 1u;
  for (unsigned tmp_size = 1; tmp_size < N; tmp_size *= 2) {
    for (unsigned i = 0, j = 0; i < 2 * tmp_size; i += 2, j++) {
      indices[i + 0] = 2 * tmp_size + 1 - tmp[j];
      indices[i + 1] = tmp[j];
    }
    memcpy(tmp, indices, sizeof(*tmp) * 2 * tmp_size);
  }

  free(tmp);
  return indices;
}

static double ComputeStep(unsigned k, unsigned N, double A, double X) {
  (void)X;
  (void)k;
  (void)A;
  // return X / N;
  return 1. / ChebyshevNode(k, N, 0., A);
}

double Scheme(unsigned N, double X, double A) {
  double h = 0.;
  const double eps = 1e-15;
  double prev = 0.;
  double current = 0.;
  double max_e = 0.;
  double solution = 0.;
  double x = 0.;
  double prod_of_norms = 1.;
  unsigned current_index = 0;
  unsigned* shuffled_indexes;

  if (fabs(log2(N) - floor(log2(N))) > eps) {
    fprintf(stderr, "failed: N should be a power of two\n");
    assert(0);
    return 0.;
  }

  shuffled_indexes = GenerateShuffledIndices(N);
  if (!shuffled_indexes) {
    fprintf(stderr, "failed: not enough memory");
    return 0.;
  }

  (void)X;
  printf("#x\tcurrent\tsolutuion\tstep\t|1-hA|\tindex\troot\n");
  for (unsigned k = 0; k <= N; k++) {
    if (k != 0) {
      current_index = shuffled_indexes[k - 1];
      h = ComputeStep(current_index, N, A, X);
    }
    current = SingleStep(k, h, A, x, prev);

    x += h;
    solution = ExactSolution(A, x);
    max_e = fmax(max_e, fabs(current - solution));
    printf("%20.15lf\t%20.15lf\t%20.15lf\t%20.15lf\t%20.15lf\t%u\t%20.15lf\n", x, current, solution, h, prod_of_norms, current_index, k == 0 ? 0. : 1. / h);
  
    if (k != N) {
      prev = current;
      prod_of_norms *= fabs(1. - h * A);
    }
  }
  printf("# (1-hA)^N: %20.15lf\n", prod_of_norms);

  free(shuffled_indexes);
  return max_e;
}
