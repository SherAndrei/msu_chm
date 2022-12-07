#pragma once

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

          double& operator[](unsigned k)       { return std::valarray<double>::operator[](k-1); }
    const double& operator[](unsigned k) const { return std::valarray<double>::operator[](k-1); }
};

double Scalar(const Vector& l, const Vector& r, double h);
