#pragma once

int Root(double *x_prev, void (*F)(double *, const double *, unsigned),
                void (*dF)(double *, const double *, unsigned), unsigned m, double eps,
                int verbose);
