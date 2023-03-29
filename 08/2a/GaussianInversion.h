#pragma once

// Inverse A to A_inversed
// A is expected to be stored by rows
int GaussMaxCol(double *A, double *A_inversed, const unsigned N);
