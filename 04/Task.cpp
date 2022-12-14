#include "Task.h"

#include <cmath>

double ExactSolution(double x)
{
    return M_PI*x + std::sin(M_PI*x);
}

double Addendum(double x)
{
    (void)x;
    return 1.;
}

double RightPart(double x)
{
    return M_PI*M_PI*std::sin(M_PI*x) + Addendum(x)*ExactSolution(x);
}
