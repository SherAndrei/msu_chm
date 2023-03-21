#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Next difference eq:
//  \frac{y_{k+1}-2y_k+y_{k-1}}{h^2} = -\lambda y_k,\ 1 \leq k \leq N-1
//  y_0 = 0
//  y_N = y_{N-1}
// for given N there are N-1 eigen vectors in N-1XN-1 matrix

static void Usage(const char *prog_name) {
  printf("Usage: %s N [m]\n"
         "\tunsigned N - matrix dim, N > 1\n"
         "\tunsigned m - num of desired eigen vector\n"
         "\t\t 0 < m < N, as the amount of eigen vectors with given initial "
         "conditions is N-1\n"
         "\tIf m is unspecified:\n"
         "\t\tPrints minimal scalar product between eigen vectors\n"
         "\telse:\n"
         "\t\t1. Prints squared length of the desired vector\n"
         "\t\t2. Substitues corresponding eigen value to the original problem "
         "and prints result\n",
         prog_name);
}

static unsigned ParseToUnsigned(const char *param_name, const char *from,
                                unsigned def) {
  unsigned ret = def;
  if (sscanf(from, "%u", &ret) != 1) {
    printf("Error parsing %s has occured, setting default value %u\n",
           param_name, def);
  }
  return ret;
}

static void FillEigenVector(unsigned m, unsigned N, double *y) {
  const double normalization_constant = sqrt(2.);
  const double fraction_in_angle = M_PI * (2. * m - 1.) / (2. * N - 1.);
  for (unsigned k = 1u; k <= N - 1; ++k) {
    y[k] = normalization_constant * sin(fraction_in_angle * k);
  }
}

static double Lambda(unsigned N, double h, unsigned m) {
  const double angle = M_PI * (2. * m - 1.) / (2 * (2. * N - 1.));
  return -4. * pow(h, -2.) * pow(sin(angle), 2.);
}

static double Scalar(const double *l, const double *r, unsigned N, double h) {
  double ret = 0.;
  for (unsigned k = 1u; k <= N - 1; k++)
    ret += l[k] * r[k];
  return ret * h;
}

static double SubstituteToTheProblem(unsigned N, const double *y, unsigned m,
                                     double h) {
  const double rev_denom = pow(h, -2.);
  const double lamda = Lambda(N, h, m);
  double left_part_numerator;
  double res;
  double *er = malloc(N * sizeof(*er));
  if (!er) {
    fprintf(stderr, "Not enough memory\n");
    return 0.;
  }

  for (unsigned k = 1u; k <= N - 1; ++k) {
    if (k == 1)
      left_part_numerator = (y[k + 1] - (2. * y[k]));
    else if (k == N - 1)
      left_part_numerator = (-y[k] + y[k - 1]);
    else
      left_part_numerator = (y[k + 1] - (2. * y[k]) + y[k - 1]);
    er[k] = left_part_numerator * rev_denom - lamda * y[k];
  }

  res = sqrt(Scalar(er, er, N, h)) / fabs(lamda);
  free(er);
  return res;
}

static double MaxScalar(unsigned N, double h) {
  double max_scalar = 0.;
  double abs_scalar;
  double *y1 = malloc(N * sizeof(*y1));
  double *y2 = malloc(N * sizeof(*y2));
  if (!y1 || !y2) {
    fprintf(stderr, "Not enough memory\n");
    free(y1);
    free(y2);
    return 0.;
  }

  for (unsigned m_i = 1u; m_i <= N - 1; ++m_i) {
    FillEigenVector(m_i, N, y1);
    for (unsigned m_j = m_i; m_j <= N - 1; ++m_j) {
      if (m_i == m_j)
        continue; // skip length
      FillEigenVector(m_j, N, y2);
      abs_scalar = fabs(Scalar(y1, y2, N, h));
      if (abs_scalar > max_scalar) {
        max_scalar = abs_scalar;
      }
    }
  }

  free(y1);
  free(y2);
  return max_scalar;
}

int main(int argc, const char *argv[]) {
  unsigned N, m = 0;
  double h;
  double *desired_vec = NULL;

  if (argc < 2 || argc > 3) {
    Usage(argv[0]);
    return 1;
  }

  N = ParseToUnsigned("N", argv[1], 2);
  if (N < 2) {
    printf("N should be greater than 1\n");
    return 2;
  }

  if (argc == 3) {
    m = ParseToUnsigned("m", argv[2], 0);
    if (!(0 < m && m < N)) {
      printf("m (%d) should be between 0 and N (%d)\n", m, N);
      return 3;
    }
  }

  h = 2. / (2. * N - 1.);

  if (m == 0) {
    printf("Maximal scalar product = %e\n", MaxScalar(N, h));
    return 0;
  }

  desired_vec = malloc(N * sizeof(*desired_vec));
  if (!desired_vec) {
    fprintf(stderr, "Not enough memory\n");
    return 4;
  }

  FillEigenVector(m, N, desired_vec);

  printf("Length squared: %e\n", Scalar(desired_vec, desired_vec, N, h));
  printf("Result of the problem with substituted eigen value: %e\n",
         SubstituteToTheProblem(N, desired_vec, m, h));
  free(desired_vec);
  return 0;
}
