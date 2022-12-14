#include "Vector.h"

#include <numeric>

Vector FillVector(double h, unsigned N, double (*filler)(double))
{
    Vector res{N};
    std::iota(begin(res), end(res), 1.);
    return (res * h).apply(filler);
}

double Scalar(const Vector& l, const Vector& r, double h)
{
    return (l * r * h).sum();
}
