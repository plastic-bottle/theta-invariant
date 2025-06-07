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

/* Puts float matrix A into REF and stores in result without scaling any rows to preserve determinant */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions */
void row_echelon_form(struct float_matrix* const A, struct float_matrix* const result) 
{
    /* Copy entries of A into result if they are different */
    if (A != result) {
        for (size_t i = 0; i < A->rows; i++) {
            for (size_t j = 0; j < A->cols; j++) {
                MATRIX_ELEMENT(result, i, j) = MATRIX_ELEMENT(A, i, j);
            }
        }
    }

    size_t current_col = 0;
    THETA_FLOAT pivot, multi;
    for (size_t i = 0; i < result->rows; i++) {
        /* Advance current_col to the next nonzero entry in the row */
        while (MATRIX_ELEMENT(result, i, current_col) == 0) {
            current_col++;
        }

        /* First element of row i */
        pivot = MATRIX_ELEMENT(result, i, current_col);
        for (size_t j = 0; j < result->rows; j++) {
            /* Skip row if its the same row as the pivot of if leading entry is 0 */
            if (j == i || MATRIX_ELEMENT(result, j, current_col) == 0) {
                /* No work to do for this row */
                continue;
            }

            /* Subtract multi * row i from row j */
            multi = MATRIX_ELEMENT(result, j, current_col) / pivot;
            for (size_t k = current_col; k < result->cols; k++) {
                MATRIX_ELEMENT(result, j, k) -= multi * MATRIX_ELEMENT(result, i, k);
            }
        }

        /* We can always do this because we have just cleared a column to contain all 0's except one row */
        current_col++;
    }
}