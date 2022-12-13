#include "Solve.h"
#include "Vector.h"

#include <cassert>

Vector Solve(const Vector& P, const Vector& F, double h, unsigned N)
{
    Vector alpha(N), beta(N), yn(N);

    const double c1 = 2./(h*h) + P[1];
    const double b1 = 1./(h*h);
    const double f1 = F[1];

    alpha[2] = b1/c1;
    beta[2] = f1/c1;

    for (auto k = 2u; k <= N-2; k++)
    {
        const double bk = b1;
        const double ck = 2./(h*h)+P[k];
        const double ak = b1;
        const double fk = F[k];

        alpha[k+1] = bk/(ck-alpha[k]*ak);
        beta[k+1]  = (alpha[k]*beta[k]+fk)/(ck-alpha[k]*ak);
    }

    assert((std::abs(alpha) < Vector{N, 1.}).min() && "Tridiagonal method is unstable!");

    const double last_c = 1./(h*h)+P[N-1];
    const double last_a = b1;

    yn[N-1] = (F[N-1]+alpha[N-1]*beta[N-1])/(last_c-last_a*alpha[N-1]);
	
    for (auto k = N - 2; k >= 1; k--)
        yn[k] = alpha[k + 1] * yn[k + 1] + beta[k + 1];

    return yn;
}
