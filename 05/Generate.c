#include <stdio.h>
#include <stdlib.h>

static void Usage(const char* argv0)
{
    printf(
        "Usage: %s N left right [output.txt]\n"
    	"\tunsigned N - amount of points to generate, N > 2\n"
        "\tdouble left - left bound of the segment\n"
        "\tdouble right - right bound of the segment\n"
		"\tfilename: output file, default -- stdout\n"
        , argv0);
}

static FILE* OpenFile(int argc, const char* filename)
{
    int file_pos = 5;
	if (argc == file_pos - 1)
		return stdout;
	return fopen(filename, "w");
}

static unsigned ParseToUnsigned(const char* param_name, const char* from) {
	unsigned ret = 0;
	if (sscanf(from, "%u", &ret) != 1) {
		fprintf(stderr, "Error parsing %s has occured, terminate\n", param_name);
        abort();
	}
	return ret;
}

static double ParseToDouble(const char* param_name, const char* from) {
	double ret = 0.;
	if (sscanf(from, "%lf", &ret) != 1) {
		fprintf(stderr, "Error parsing %s has occured, terminate\n", param_name);
        abort();
	}
	return ret;
}

static void GenerateX(unsigned N, double left, double right, double* x)
{
    double step = (right - left) / (N-1.);
    x[0] = left;
    for (int i = 0; i < N; i++)
    {
        x[i+1] = x[i] + step;
    }
}

static double GenerateYCoordinate(double x)
{
    return x * x;
}

static void GenerateY(const double* x, unsigned N, double* y)
{
    for (int i = 0; i < N; i++)
        y[i] = GenerateYCoordinate(x[i]);
}

int main(int argc, const char* argv[])
{
    FILE* out = NULL;
    double left_bound = 0.;
    double right_bound = 0.;
    double* x = NULL;
    double* y = NULL;
    unsigned N = 0;

    if (argc < 4 || argc > 5) {
        Usage(argv[0]);
        return 1;
    }

	N = ParseToUnsigned("N", argv[1]);
    left_bound = ParseToDouble("left bound", argv[2]);
    right_bound = ParseToDouble("right bound", argv[3]);

    out = OpenFile(argc, argv[4]);
	if (!out) {
		fprintf(stderr, "Cannot open output file\n");
		return 3;
	}

    x = (double*)malloc(sizeof(double)*N);
    y = (double*)malloc(sizeof(double)*N);
    if (!x || !y)
    {
        fprintf(stderr, "Not enough memory\n");
        fclose(out);
        free(x);
        free(y);
        return 4;
    }  

    GenerateX(N, left_bound, right_bound, x);
    GenerateY(x, N, y);

    fprintf(out, "%u\n", N);
    for (int i = 0; i < N; i++)
    {
        fprintf(out, "%lf %lf\n", x[i], y[i]);
    }

    fclose(out);
    free(x);
    free(y);
}
