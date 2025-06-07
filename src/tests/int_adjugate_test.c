/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of theta_invariant.                                     *
 *                                                                            *
 *  theta_invariant is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  theta_invariant is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License along   *
 *  with theta_invariant. If not, see <https://www.gnu.org/licenses/>.        *
 ******************************************************************************/

#include "../theta_implementation.h"

void int_adjugate_test() {
    int n = 19; // odd
    THETA_FLOAT t = (n + 1) / 2.0;
    int N = 2 * n + 1;

    struct float_matrix* A = make_float_matrix(N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            MATRIX_ELEMENT(A, i, j) = 0;
        }
    }
    for (int i = 0; i < n; i++) {
        if (i % 2 == 0) {
            MATRIX_ELEMENT(A, i, i+1) = -t;
            MATRIX_ELEMENT(A, i, n+i+1) = t-1;
            MATRIX_ELEMENT(A, n+i, n+i+1) = -1;
        } else {
            MATRIX_ELEMENT(A, n+i, n+i+1) = -t;
            MATRIX_ELEMENT(A, n+i, i+1) = t-1;
            MATRIX_ELEMENT(A, i, i+1) = -1;
        }
    }
    for (int i = 0; i < N; i++) {
        MATRIX_ELEMENT(A, i, i) += 1;
    }

    struct int_matrix* result = make_int_matrix(N, N);

    THETA_INT determinant = int_adjugate(A, result);

    /*for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f   ", (double) MATRIX_ELEMENT(A, i, j));
        }
        printf("\n");
    }*/

    printf("\n%I64d\n\n", determinant);

    /*for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f\n", (double) MATRIX_ELEMENT(result, i, j));
        }
        printf("\n");
    }*/
}