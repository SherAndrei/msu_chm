#include "Exact.h"
#include "Vector.h"

#include <numeric>

Vector Exact(double, unsigned N)
{
    return Vector{N};
}

double EigenValue(unsigned m, double h, unsigned N)
{
    const double angle = M_PI * (2. * m - 1.) / (2 * (2. * N - 1.));
    return std::pow(2. * std::sin(angle) / h, 2.);
}

Vector EigenVector(unsigned m, unsigned N)
{
    Vector y(N);
    std::iota(begin(y), end(y), 1.);
    const double normalization_constant = std::sqrt(2.);
    const double fraction_in_angle = M_PI * (2. * m - 1.) / (2. * N - 1.);
    return { normalization_constant * std::sin(y * fraction_in_angle) };
}
