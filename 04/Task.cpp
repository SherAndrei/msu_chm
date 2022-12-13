#include "Task.h"
#include "Vector.h"

#include <numeric>

static Vector FillSegment(double step, unsigned amount, double(*filler)(double))
{
    Vector res{ amount };
    std::iota(begin(res), end(res), 1.);
    return Vector{ (res * step).apply(filler) };
}

static double ExactSolutionValue(double x)
{
    return x*x - 2*x;
}

Vector ExactSolution(double h, unsigned N)
{
    return FillSegment(h, N, &ExactSolutionValue);
}

static double AddendumValue(double x)
{
    (void)x;
    return 1.;
}

Vector Addendum(double h, unsigned N)
{
    return FillSegment(h, N, &AddendumValue);
}

static double RightPartValue(double x)
{
    return -2 + ExactSolutionValue(x);
}

Vector RightPart(double h, unsigned N)
{
    return FillSegment(h, N, &RightPartValue);
}
