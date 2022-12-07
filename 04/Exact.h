#pragma once

struct Vector;

Vector Exact(double h, unsigned N);

double EigenValue(unsigned m, double h, unsigned N);
Vector EigenVector(unsigned m, unsigned N);
