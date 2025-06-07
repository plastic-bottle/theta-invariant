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

#include "theta_implementation.h"

/* Test */
int main() {
    struct float_matrix* A = (struct float_matrix*) safe_malloc(sizeof(struct float_matrix));
    A->rows = 3;
    A->cols = 3;
    A->data = (THETA_FLOAT*) safe_malloc(A->rows * A->cols * sizeof(THETA_FLOAT));

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            MATRIX_ELEMENT(A, i, j) = (i+2) * j - 1.0;
        }
    }

    struct float_matrix* result = (struct float_matrix*) safe_malloc(sizeof(struct float_matrix));
    result->rows = 3;
    result->cols = 3;
    result->data = (THETA_FLOAT*) safe_malloc(A->rows * A->cols * sizeof(THETA_FLOAT));

    row_echelon_form(A, result);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f   ", MATRIX_ELEMENT(A, i, j));
        }
        printf("\n");
    }

    return 0;
}