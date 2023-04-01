#pragma once

void InitialApproximation(double *x, unsigned m);

void F(double *f, const double *x, unsigned m);

void dF(double *jacobian, const double *x, unsigned m);
