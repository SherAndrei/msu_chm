
#include "Matrix.h"

#include <stdio.h>
#include <stdlib.h>

static void Usage(const char* argv0)
{
    printf(
        "Usage: %s [input]\n"
        "\tCalculate interpolation polynom by using input data.\n"
        "\tResult is printed to output.txt\n"
        "Parameters:\n"
    	"\tunsigned N - amount of points to generate, N > 2\n"
        "\tdouble left - left bound of the segment\n"
        "\tdouble right - right bound of the segment\n"
        "\tinput: input file, default -- stdin\n"
        "\t\tExpected format:\n"
        "\t\tN\n"
        "\t\tx1\ty1\n"
        "\t\t...\t...\n"
        "\t\tx_{N-1}\ty_{N-1}\n"
		"\toutput: output file, default -- stdout\n"
        "\t\tOutput format:\n"
        "\t\tX Exact Gaussian Delta\n"
        , argv0);
}

static FILE* OpenInputFile(int argc, const char* argv[])
{
    int file_pos = 1;
	if (argc == file_pos)
		return stdin;
	return fopen(argv[1], "r");
}

static void SolveWithGaussianElimination(const double* x, const double* y, unsigned N, double* a)
{
    const double eps = 1e-15;
    double norm = 0.;
    struct Matrix A;
    int res = 0;
    InitMatrix(x, N, &A);

    norm = MatrixNorm(A);
    res = GaussMaxCol(y, norm * eps, A, a);

    if (res == -1) {
        fprintf(stderr, "Matrix is degenerate!\n");
    }
    DestroyMatrix(&A);
}

int main(int argc, const char* argv[])
{
    FILE* in = NULL;
    // FILE* out = NULL;
    int error = 0;
    unsigned N = 0;
    double* x = NULL;
    double* y = NULL;
    double* a = NULL;

    if (argc > 2) {
        Usage(argv[0]);
        return 1;
    }
    
    in = OpenInputFile(argc, argv);
    if (!in) {
        fprintf(stderr, "Cannot open input file\n");
        return 2;
    }

    if (fscanf(in, "%u", &N) != 1) {
        fprintf(stderr, "Error parsing N from input\n");
        error = 3;
        goto clear;
    }

    x = (double*)malloc(N*sizeof(double));
    y = (double*)malloc(N*sizeof(double));
    a = (double*)malloc(N*sizeof(double));
    if (!x || !y || !a) {
        fprintf(stderr, "Not enough memory\n");
        error = 4;
        goto clear;
    }

    for (unsigned i = 0u; i < N; i++) {
        if (fscanf(in, "%lf%lf", x + i, y + i) != 2) {
            fprintf(stderr, "Error parsing input data\n");
            error = 4;
            goto clear;
        }
    }

    SolveWithGaussianElimination(x, y, N, a);

clear:
    fclose(in);
    free(x);
    free(y);
    free(a);
    return error;
}
