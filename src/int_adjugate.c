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

/* Returns the determinant of A, and populates matrix result with the entries of the adjugate of A */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions and are both square */
THETA_INT int_adjugate(struct int_matrix* const A, struct int_matrix* const result)
{
    size_t N = A->rows;

    /* Set up the n by 2n augmented matrix that will be used to compute inverse */
    struct float_matrix* augment = make_float_matrix(N, 2 * N);

    /* Copy elements of A into augment to set it up */
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            MATRIX_ELEMENT(augment, i, j) = (THETA_FLOAT) MATRIX_ELEMENT(A, i, j);
        }
    }

    /* Set up the right half of augment with the identity matrix */
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            MATRIX_ELEMENT(augment, i, j + N) = ((i == j) ? 1 : 0);
        }
    }

    /* Run REF on augment */
    row_echelon_form(augment, augment);

    /* Compute determinant */
    THETA_FLOAT temp_determinant = 1;
    for (size_t i = 0; i < N; i++) {
        temp_determinant *= MATRIX_ELEMENT(augment, i, i);
    }
    THETA_INT determinant = (THETA_INT) (temp_determinant + 0.5);

    /* Divide each row by the diagonal entry in the left half and then multiply by the determinant */
    /* This will result in the identity matrix on the left and the adjugate matrix on the right */
    THETA_FLOAT scale;
    for (size_t i = 0; i < N; i++) {
        /* Number to divide by */
        scale = MATRIX_ELEMENT(augment, i, i);
        if (scale == 0) {
            exit(EXIT_FAILURE);
        }

        /* Loop through columns to divide by scale and multiply by determinant */
        for (size_t j = 0; j < N; j++) {
            MATRIX_ELEMENT(augment, i, j + N) = MATRIX_ELEMENT(augment, i, j + N) / scale * determinant;
        }
    }

    /* Write adjugate to result */
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            MATRIX_ELEMENT(result, i, j) = (THETA_INT) (MATRIX_ELEMENT(augment, i, j + N) + 0.5);
        }
    }

    safe_free(augment->data);
    safe_free(augment);

    return determinant;
}
