#include "Solve.h"
#include "Task.h"
#include "Vector.h"

#include <numeric>

static double EigenValue(const Vector& p, unsigned m, double h, unsigned N)
{
    const double angle = M_PI * (2. * m - 1.) / (2 * (2. * N - 1.));
    return std::pow(2. * std::sin(angle) / h, 2.) + p[m];
}

static Vector EigenVector(unsigned m, unsigned N)
{
    Vector y(N);
    std::iota(begin(y), end(y), 1.);
    const double normalization_constant = std::sqrt(2.);
    const double fraction_in_angle = M_PI * (2. * m - 1.) / (2. * N - 1.);
    return { normalization_constant * std::sin(y * fraction_in_angle) };
}

Vector Solve(double h, unsigned N)
{
    Vector y(N);
    auto f = FillVector(h, N, &RightPart);
    auto p = FillVector(h, N, &Addendum);
    for (auto m = 1u; m <= N-1; m++) 
    {
        const Vector em = EigenVector(m, N);
        const double dm = Scalar(f, em, h);
        const double cm = dm / EigenValue(p, m, h, N);
        y += cm * em;
    }
    return y;
}
