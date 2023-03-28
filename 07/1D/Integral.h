#pragma once

double IntegralOnSegment(double a, double b, double (*f)(double));

inline double Integral(double a, double b, double (*f)(double), unsigned N) {
  const double step = (b - a) / N;
  double sum = 0.;
  for (unsigned i = 0; i < N; ++i) {
    sum += IntegralOnSegment(a, a + step, f);
    a += step;
  }
  return sum;
}
