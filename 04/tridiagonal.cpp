#include "Solve.h"
#include "Vector.h"

Vector Solve(double h, unsigned N)
{
    (void)h;
    return Vector{N};
#if 0
    Vector alpha(N), beta(N), res(N);
	const double H_2=h*h;
    
    // alpha[1] = (2./H_2)/(2./H_2+ b(0));
    // beta[1] = f(0)/(2./H_2 + b(0));
    
    // for (int k = 1; k <= N - 1; k++) 
	// {
    //     alpha[k + 1] = (1./H_2)/(2./H_2 + b(k*h) - alpha[k]/H_2);
    //     beta[k + 1] = (f(k*h) + beta[k]/H_2)/(2./H_2 + b(k*h) - alpha[k]/H_2);
    // }
    // y[N] = (f(N*h) + 2*beta[N]/H_2)/(2./H_2 + b(N*h) - 2*alpha[N]/H_2);
	
    // for (int k = N - 1; k >= 0; k--) y[k] = alpha[k + 1]*y[k + 1] + beta[k + 1];
  
    return res;
#endif
}
