#include "step.h"

#include <cmath>
#include <cstdio>

void E_n(int m, double A);
double ExactSolution(double A, double h, double k);
void Usage(const char* argv0);

double ExactSolution(double A, double x_k)
{
	return std::exp(-A * x_k);
}

void E_n(int m, double A)
{
	double h = 0.1;
	double el;
	double next;
    double max_e = 0.;
    double razn = 0.;
	int n;
    int N = 10;
	int fl=0;

    for(n = 1; n <= 3; n++)
	{
        for(int k = 0; k<=N; k++)
		{
/*            if(h*A >= 1.)
			{
                printf("E%d = Err             ", n);
				fl=1;
                break;
			}*/

            Step(k, h, A, el, next);

            if(fabs(next) > 1e200)
		    {
                max_e = -1;
                break;
			}
            if(fabs(next) < 1e-100)
			{
                max_e = 0;
                break;
			}
            razn = next - ExactSolution(A, (k * h));
            razn = fabs(razn);
            if(razn > max_e) max_e = razn;
		}
        if(fl==0)printf ("E%d = %e    ", n, max_e);
        h /= 10;
        N *= 10;
		fl=0;
    }

	h /= 100;
	N *= 100;
	n = 6;

	for(int k = 0; k<=N; k++)
	{
		/*if(h*A >= 1.)
		{
			printf("E%d = Err", n);
			break;
		}*/

		Step(k, h, A, el, next);
		if(fabs(next) > 1e200)
		{
			max_e = -1.;
			break;
		}
		if(fabs(next) < 1e-100)
		{
			max_e = 0.;
			break;
		}
		razn = next - exp(-A*(k*h));
		razn = fabs(razn);
		if (razn > max_e) max_e = razn;
	}
	printf("E%d = %e    ", n, max_e);
	h /= 10;
	N *= 10;
	
    printf("m = %d    A = %lf, h=%e\n", m, A, h);
}

void Usage(const char* argv0) {
    std::printf(
        "Usage: %s m A\n"
		"\tunsigned m - order of convergence\n"
        "\tdouble A - task parameter\n"
        , argv0);
}

int main(int argc, const char* argv[])
{
	if (argc != 3) {
		Usage(argv[0]);
		return 1;
	}

	double A = 0.;
	unsigned m = 0;
	if (!(std::sscanf(argv[1], "%lf", &A) == 1
	   && std::sscanf(argv[1], "%u", &m) == 1)) {
		std::printf("Error: cannot parse input parameters\n");
		return 2;
	}

	E_n(m, A);
}
