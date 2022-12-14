#pragma once

#include <cassert>
#include <valarray>

struct Vector : public std::valarray<double>
{
public:
    explicit Vector(unsigned N, double def = 0.)
        : std::valarray<double>(def, N-1) {}

    Vector(const std::valarray<double>& other)
        : std::valarray<double>(other) {}
    Vector(std::valarray<double>&& other)
        : std::valarray<double>(std::move(other)) {}

    double& operator[](unsigned k) {
        assert_on_incorrect_index(k);
        return std::valarray<double>::operator[](k-1);
    }
    const double& operator[](unsigned k) const {
        assert_on_incorrect_index(k);
        return std::valarray<double>::operator[](k-1);
    }

private:
    void assert_on_incorrect_index([[maybe_unused]] unsigned k) const {
        [[maybe_unused]] const auto N = size() + 1;
        assert(1 <= k && k <= N-1);
    }
};

Vector FillVector(double h, unsigned N, double (*filler)(double));
double Scalar(const Vector& l, const Vector& r, double h);
