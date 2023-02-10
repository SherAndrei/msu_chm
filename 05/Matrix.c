#include "Matrix.h"

#include <math.h>
#include <stdlib.h>

#define A(i, j) (A[(i) * N + (j)])

int InitMatrix(const double* x, const unsigned N, struct Matrix* m)
{
    double* const A = m->A = (double*)malloc(N*N*sizeof(double));
    if (!A)
        return -1;

    for (unsigned j = 0u; j < N; j++)
        for (unsigned i = 0u; i < N; i++)
            A(j, i) = pow(x[i], j);
    m->N = N;
    return 0;
}

void PrintMatrix(struct Matrix m) {
    const double* const A = m.A;
    const unsigned N = m.N;
    for (unsigned i = 0u; i < N; ++i) {
        for (unsigned j = 0u; j < N; ++j)
            printf("%e%c", A(j,i), " \n"[j == N - 1]);
    }
}

void DestroyMatrix(struct Matrix* m)
{
    free(m->A);
    m->N = 0;
}

double MatrixNorm(const struct Matrix m)
{
    const double* const A = m.A;
    const unsigned N = m.N;
    unsigned i, j;
    double max = 0., current = 0.;
    for (i = 0; i < N; i++) {
        current = 0.;
        for (j = 0; j < N; j++) {
            current += fabs(A(j, i));
        }
        if (current > max)
            max = current;
    }
    return max;
}

static inline void Swap(double* const lhs, double* const rhs) {
    double temp = *lhs;
    *lhs = *rhs;
    *rhs = temp;
}

// Solve Aa = y
int GaussMaxCol(const double* y, double error, const struct Matrix m, double* a)
{
    double* const A = m.A;
    const unsigned N = m.N;
    unsigned i, j, k;
    double max_elem_in_columns = 0.;
    unsigned column_with_max_elem = 0;
    double c = 0.;

    for (i = 0; i < m.N; i++) {
        a[i] = y[i];
    }

    for (j = 0; j < m.N; j++) {
        max_elem_in_columns = fabs(A(j, j));
        column_with_max_elem = j;

        // find max elem between columns
        for (i = j; i < m.N; i++) {
            if (fabs(A(i, j)) > max_elem_in_columns) {
                column_with_max_elem = i;
                max_elem_in_columns = fabs(A(i, j));
            }
        }

        if (max_elem_in_columns < error) {
            // matrix is degenerate!
            return -1;
        }

        // swap current row with row with max elem
        if (column_with_max_elem != j) {
            for (i = j; i < m.N; i++) {
                Swap(&(A(j, i)), &(A(column_with_max_elem, i)));
            }
            Swap(a + j, a + column_with_max_elem);
        }

        // divide all row by first elem
        c = 1. / A(j, j);
        A(j, j) = 1.;
        for (i = j + 1; i < m.N; i++) {
            A(j, i) *= c;
        }
        a[j] *= c;

        for (i = j + 1; i < m.N; i++) {
            c = A(i, j);
            if (fabs(c) > error) {
                A(i, j) = 0.;
                for (k = j + 1; k < m.N; k++) {
                    A(i, k) -= c * A(j, k);
                }
            }
        }
        PrintMatrix(m);
        printf("\n\n");
    }
    // reverse step
    // for (j = m.N - 1; j >= 0; j--) {
    //     for (i = 0; i < j; i++) {
    //         c = A(i, j);
    //         if (fabs(c) > error) {
    //             for (m = 0; m < n; m++)
    //                 E(i, m) -= c * E(j, k);
    //             A(i, j) = 0.;
    //         }
    //     }
    // }
    return 0;
}
