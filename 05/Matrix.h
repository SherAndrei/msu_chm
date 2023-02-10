#pragma once

#include <stdio.h>

struct Matrix {
    double* A;
    unsigned N;
};

int InitMatrix(const double* x, const unsigned N, struct Matrix* m);
void PrintMatrix(struct Matrix m);
void DestroyMatrix(struct Matrix* m);

double MatrixNorm(const struct Matrix m);
int GaussMaxCol(const double* y, double error, const struct Matrix m, double* a);
