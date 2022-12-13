#include "Task.h"
#include "Vector.h"

#include <numeric>

static double ExactSolutionValue(double x)
{
    return x*x - 2*x;
}

Vector ExactSolution(double h, unsigned N)
{
    Vector res{N};
    std::iota(begin(res), end(res), 1.);
    return Vector{ (res * h).apply(&ExactSolutionValue) };
}

Vector Addendum(double h, unsigned N)
{
    (void)h;
    return Vector{ N, 1. };
}

static double RightPartValue(double x)
{
    return x*x-2*x-2;
}

Vector RightPart(double h, unsigned N)
{
    Vector res{N};
    std::iota(begin(res), end(res), 1.);
    return Vector{ (res * h).apply(&RightPartValue) };
}
