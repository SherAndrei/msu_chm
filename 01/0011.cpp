#include <cmath>
#include <cstdio>
#include <vector>

// Next difference eq:
//  \frac{y_{k+1}-2y_k+y_{k-1}}{h^2} = -\lambda y_k,\ 1 \leq k \leq N-1
//  y_0 = 0
//  y_N = y_{N-1}
// for given N there are N-1 eigen vectors in N-1XN-1 matrix

class Vector
{
public:
    Vector(unsigned N)
        : m_elems(N-1) {}

          double& operator[](unsigned k)       { return m_elems[k-1]; }
    const double& operator[](unsigned k) const { return m_elems[k-1]; }

    size_t size() const { return m_elems.size(); }

private:
    std::vector<double> m_elems;
};

void Usage(const char* prog_name);
unsigned ParseToUnsigned(const char* param_name, const char* from, unsigned def);

void FillEigenVector(unsigned m, Vector& ret);
double Lambda(unsigned N, double h, unsigned m);

double Scalar(const Vector& l, const Vector& r, double h);

double SubstituteToTheProblem(unsigned N, const Vector& l, unsigned m, double h);
double MaxScalar(unsigned N, double h);

int main(int argc, const char* argv[]) {
    if (argc < 2 || argc > 3) {
        Usage(argv[0]);
        return 1;
    }
    unsigned N, m = 0;
    N = ParseToUnsigned("N", argv[1], 2);
    if (N < 2) {
        std::printf("N should be greater than 1\n");
        return 2;
    }
    if (argc == 3) {
        m = ParseToUnsigned("m", argv[2], 0);
        if (!(0 < m && m < N)) {
            std::printf("m (%d) should be between 0 and N (%d)\n", m, N);
            return 3;
        }
    }

    const double h = 2./(2.*N - 1.);

    if (m == 0) {
        std::printf("Maximal scalar product = %e\n", MaxScalar(N, h));
        return 0;
    }

    Vector desired_vec(N);
    FillEigenVector(m, desired_vec);

    std::printf("Length squared: %e\n", Scalar(desired_vec, desired_vec, h));
    std::printf("Result of the problem with substituted eigen value: %e\n", SubstituteToTheProblem(N, desired_vec, m, h));
    return 0;
}

void Usage(const char* prog_name) {
    std::printf(
        "Usage: %s N [m]\n"
        "\tunsigned N - matrix dim, N > 1\n"
        "\tunsigned m - num of desired eigen vector\n"
        "\t\t 0 < m < N, as the amount of eigen vectors with given initial conditions is N-1\n"
        "\tIf m is unspecified:\n"
        "\t\tPrints minimal scalar product between eigen vectors\n"
        "\telse:\n"
        "\t\t1. Prints squared length of the desired vector\n"
        "\t\t2. Substitues corresponding eigen value to the original problem and prints result\n"
        , prog_name);
}

unsigned ParseToUnsigned(const char* param_name, const char* from, unsigned def) {
    unsigned ret = def;
    if (std::sscanf(from, "%u", &ret) != 1) {
        std::printf("Error parsing %s has occured, setting default value %u\n", param_name, def);
    }
    return ret;
}

void FillEigenVector(unsigned m, Vector& y) {
    const auto N = y.size() + 1;
    const double normalization_constant = std::sqrt(2.);
    const double fraction_in_angle = M_PI*(2.*m-1.)/(2.*N-1.);
    for (auto k = 1u; k <= N-1; ++k) {
        y[k] = normalization_constant * std::sin(fraction_in_angle * k);
    }
}

double Lambda(unsigned N, double h, unsigned m) {
    const double angle = M_PI*(2.*m - 1.)/(2*(2.*N-1.));
    return -4. * std::pow(h, -2.) * std::pow(std::sin(angle), 2.);
}

double Scalar(const Vector& l, const Vector& r, double h) {
    double ret = 0.;
    for (auto k = 1u; k <= r.size(); k++)
        ret += l[k] * r[k] * h;
    return ret;
}

double SubstituteToTheProblem(unsigned N, const Vector& y, unsigned m, double h) {
    const auto rev_denom = std::pow(h, -2.);
    const auto lamda = Lambda(N, h, m);
    Vector er(N);
    for (auto k = 1u; k <= er.size(); ++k) {
        double left_part_numerator = 0.;
        if (k == 1)
            left_part_numerator = (y[k+1]-(2.*y[k]));
        else if (k == N - 1)
            left_part_numerator = (-y[k]+y[k-1]);
        else
            left_part_numerator = (y[k+1]-(2.*y[k])+y[k-1]);
        er[k] = left_part_numerator * rev_denom - lamda * y[k];
    }
    return std::sqrt(Scalar(er, er, h)) / std::abs(lamda);
}

double MaxScalar(unsigned N, double h) {
    double max_scalar = 0.;
    Vector y1(N), y2(N);
    for (auto m_i = 1u; m_i <= N-1; ++m_i) {
        FillEigenVector(m_i, y1);
        for (auto m_j = m_i; m_j <= N-1; ++m_j) {
            if (m_i == m_j)
                continue; // skip length
            FillEigenVector(m_j, y2);
            double abs_scalar = std::abs(Scalar(y1, y2, h));
            if (abs_scalar > max_scalar) {
                max_scalar = abs_scalar;
            }
        }
    }
    return max_scalar;
}
