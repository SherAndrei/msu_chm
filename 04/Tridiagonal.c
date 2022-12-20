#include "Solve.h"

#include <stdio.h>
#include <stdlib.h>

void Solve(const double* P, const double* F, double* Yn, double h, unsigned N)
{
    double* const alpha = (double*)calloc(N, sizeof(double));
    double* const beta = (double*)calloc(N, sizeof(double));
    const double c1 = 2./(h*h) + P[1];
    const double b1 = 1./(h*h);
    const double f1 = F[1];
    double bk;
    double ck;
    double ak;
    double fk;
    const double last_c = 1./(h*h)+P[N-1];
    const double last_a = b1;

    if (!alpha || !beta) {
        fprintf(stderr, "Not enought memory for algorithm\n");
        free(alpha);
        free(beta);
        return;
    }

    alpha[2] = b1/c1;
    beta[2] = f1/c1;

    for (unsigned k = 2u; k <= N-2; k++)
    {
        bk = b1;
        ck = 2./(h*h)+P[k];
        ak = b1;
        fk = F[k];

        alpha[k+1] = bk/(ck-alpha[k]*ak);
        beta[k+1]  = (ak*beta[k]+fk)/(ck-alpha[k]*ak);
    }

    Yn[N-1] = (F[N-1]+last_a*beta[N-1])/(last_c-last_a*alpha[N-1]);
	
    for (unsigned k = N - 2; k >= 1; k--)
        Yn[k] = alpha[k + 1] * Yn[k + 1] + beta[k + 1];

    free(alpha);
    free(beta);
}
