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

struct polynomial_matrix* multiply_polynomial_matrices(struct polynomial_matrix* A, struct polynomial_matrix* B)
{
    struct polynomial_matrix* result = make_polynomial_matrix(A->rows, B->cols);

    /* Also equal to B->rows */
    int m = A->cols;

    int* zero = {0};
    
    for (int r = 0; r < result->rows; r++) {
        for (int c = 0; c < result->cols; c++) {
            MATRIX_ELEMENT(result, r, c) = make_polynomial(0, zero);
            for (int i = 0; i < m; i++) {
                MATRIX_ELEMENT(result, r, c) = add_polynomials(
                    MATRIX_ELEMENT(result, r, c),
                    multiply_polynomials(MATRIX_ELEMENT(A, r, i), MATRIX_ELEMENT(result, i, c))
                );
            }
        }
    }

    return result;

}
