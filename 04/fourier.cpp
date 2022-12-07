#include "Exact.h"
#include "Solve.h"
#include "Vector.h"

Vector Solve(double h, unsigned N)
{
    Vector f = Exact(h, N), res(N, 0.);
    for (auto m = 0u; m <= N - 1; m++) 
    {
        Vector ym = EigenVector(m, N);
        double Cm = Scalar(ym, f, h) / EigenValue(m, h, N);
        res += Cm * ym;
    }
    return res;
}
