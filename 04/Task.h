#pragma once

struct Vector;

// -ExactSolution'' + Addendum*ExactSolution = RightPart
Vector ExactSolution(double h, unsigned N);
Vector Addendum(double h, unsigned N);
Vector RightPart(double h, unsigned N);
