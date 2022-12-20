#pragma once

// -ExactSolution'' + Addendum*ExactSolution = RightPart
void ExactSolution(double h, double* y, unsigned N);
void Addendum(double h, double* p, unsigned N);
void RightPart(double h, double* f, unsigned N);
