#include "Exact.h"
#include "Solve.h"
#include "Vector.h"

#include <cstdio>

// Next difference eq:
//  \frac{y_{k+1}-2y_k+y_{k-1}}{h^2} = -\lambda y_k,\ 1 \leq k \leq N-1
//  y_0 = 0
//  y_N = y_{N-1}
// for given N there are N-1 eigen vectors in N-1XN-1 matrix

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
    if (N == 0)
        return 2;
    
    std::FILE* out = OpenFile(argc, argv);
	if (!out) {
		fprintf(stderr, "Cannot open output file\n");
		return 3;
	}

    const double h = 2./(2.* N - 1.);
    Vector exact = Exact(h, N);
    Vector solved = Solve(h, N);

    for (auto i = 1u; i <= N; i++) {
        std::fprintf(out, "%e %e\n", solved[i], exact[i]);
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
        "\tunsigned N - matrix dim, N > 1\n"
		"\tfilename: output file, default -- stdout\n"
        "Result output is the table in format:\n"
		" yn | y \n"
        , prog_name);
}

