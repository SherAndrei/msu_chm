#include "Task.h"
#include "Solve.h"
#include "Vector.h"

#include <cstdio>

std::FILE* OpenFile(int argc, const char* argv[]);
unsigned ParseToUnsigned(const char* param_name, const char* from, unsigned def);
void Usage(const char* prog_name);

int main(int argc, const char* argv[])
{
    if (argc < 2 || argc > 3) {
        Usage(argv[0]);
        return 1;
    }
    
    unsigned N = ParseToUnsigned("N", argv[1], 0);
    if (N < 3)
        return 2;
    
    std::FILE* out = OpenFile(argc, argv);
	if (!out) {
		fprintf(stderr, "Cannot open output file\n");
		return 3;
	}

    const double h = 2./(2.* N - 1.);
    
    Vector right = RightPart(h, N);
    Vector add = Addendum(h, N);
    Vector solved = Solve(add, right, h, N);
    Vector exact = ExactSolution(h, N);

    std::fprintf(out, "%e\t%e\t%e\n", 0., 0., 0.); // initial condition
    for (auto i = 1u; i <= N-1; i++) {
        std::fprintf(out, "%e\t%e\t%e\n", i*h, solved[i], exact[i]);
    }

    std::printf("Error: %e\n", std::sqrt(std::pow((exact - solved), Vector{N, 2.}).sum() / N));

    fclose(out);
}

std::FILE* OpenFile(int argc, const char* argv[])
{
	if (argc == 2)
		return stdout;
	return std::fopen(argv[2], "w");
}

unsigned ParseToUnsigned(const char* param_name, const char* from, unsigned def) {
    unsigned ret = def;
    if (std::sscanf(from, "%u", &ret) != 1) {
        std::fprintf(stderr, "Error parsing %s has occured, setting default value %u\n", param_name, def);
    }
    return ret;
}

void Usage(const char* prog_name) {
    std::printf(
        "Usage: %s N [output.txt]\n"
        "\tunsigned N - matrix dim, N > 2\n"
		"\tfilename: output file, default -- stdout\n"
        "Solving difference equation -y''(x)+P(x)y(x)=F(x)\n"
        "with initial conditions y(0)=0, y'(1)=0 using next\n"
        "difference scheme -\\frac{y_{k-1}-2y_k+y_{k+1}}{h^2}+P_ky_k=F_k\n"
        "and corresponding initial conditions y_0=0, y_N=y_{N-1}.\n"
        "Result output is the table in format:\n"
		"x | yn | y \n"
        , prog_name);
}
