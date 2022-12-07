#include "Vector.h"

double Scalar(const Vector& l, const Vector& r, double h)
{
    return (l * r * h).sum();
}
