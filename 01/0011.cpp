#include <cassert>
#include <cmath>
#include <cstdio>
#include <charconv>
#include <vector>

void Usage(const char* prog_name);
double EigenVectorCoordinate(unsigned k, unsigned N, unsigned m);
double EigenValue(unsigned N, double h, unsigned m);

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        Usage(argv[0]);
        return 1;
    }
    unsigned N;
    if (std::sscanf(argv[1], "%u", &N)) {
        std::puts("Invalid N");
        return 2;
    }
    std::vector<double> matr(N);
    const double h = 1./N;
    (void)h;
}

void Usage(const char* prog_name) {
    std::printf("Usage: %s <N>\n", prog_name);
}

double EigenVectorCoordinate(unsigned k, unsigned N, unsigned m) {
    using std::pow;
    assert(k <= N);
    assert(m >= 1 && m <= N-1);
    const double angle = M_PI*(2.*m - 1.)*k/(2.*N-1.);
    return std::sin(angle);
}

double EigenValue(unsigned N, double h, unsigned m) {
    using std::pow;
    assert(m >= 1 && m < N - 1);
    const double angle = M_PI*(2.*m - 1.)/(2*(2.*N-1.));
    return 4. * pow(h, -2.) * pow(std::sin(angle), 2.);
}

