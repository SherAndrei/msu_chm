#pragma once

unsigned NumberOfEquations();

// assumes y = malloc(sizeof(double) * NumberOfEquations())
void ExactSolution(double* y, double x);

// assumes:
//  y = malloc(sizeof(double) * NumberOfEquations())
//  f = malloc(sizeof(double) * NumberOfEquations())
void RightPartOfEquations(double* f, const double* y, double x);
