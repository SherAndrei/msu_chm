#include "Solve.h"
#include "Task.h"
#include "Vector.h"

#include <cassert>

Vector Solve(double h, unsigned N)
{
    Vector alpha(N), beta(N), yn(N);

    const double h_squared = h*h;
    const double c1 = 2./h_squared + Addendum(h);
    const double b1 = 1./h_squared;
    const double f1 = RightPart(h);

    alpha[2] = b1/c1;
    beta[2] = f1/c1;

    for (auto k = 2u; k <= N-2; k++)
    {
        const double bk = b1;
        const double ck = 2./h_squared+Addendum(h*k);
        const double ak = b1;
        const double fk = RightPart(h*k);

        alpha[k+1] = bk/(ck-alpha[k]*ak);
        beta[k+1]  = (alpha[k]*beta[k]+fk)/(ck-alpha[k]*ak);
    }

    assert((std::abs(alpha) < Vector{N, 1.}).min() && "Tridiagonal method is unstable!");

    const double last_c = 1./h_squared+Addendum(h*(N-1));
    const double last_a = b1;

    yn[N-1] = (RightPart(h*(N-1))+alpha[N-1]*beta[N-1])/(last_c-last_a*alpha[N-1]);
	
    for (auto k = N - 2; k >= 1; k--)
        yn[k] = alpha[k + 1] * yn[k + 1] + beta[k + 1];

    return yn;
}
