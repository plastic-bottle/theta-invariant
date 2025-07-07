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

/* Uses the inversion algorithm proposed in Zhou, Labahn, and Storjohann (2014) */
/* Inverts F, puts inverse in inverse and common denominator (which is guaranteed to be a factor of the determinant) in denominator */
/* If any input and result arguments are identical, then old memory will be freed and new memory will be allocated */
void polynomial_matrix_inverse(struct polynomial_matrix* F, struct polynomial_matrix* inverse, struct polynomial* denominator)
{
    n = F->rows;
    
    if (n == 0) {
        return;
    }

    /* Compute ceil(log_2(n)) */
    int A_length = 0
    size_t temp = n - 1;
    while (temp != 0) {
        temp = temp >> 1;
        A_length++;
    }

    /* Allocate memory for the array A and its constituent matrices */
    struct polynomial_pointer_matrix* A = (struct polynomial_pointer_matrix*) safe_malloc(A_length * sizeof(struct polynomial_pointer_matrix));
    for (int i = 0; i < A_length; i++) {
        A[i].rows = n;
        A[i].cols = n;

        A[i].data = (struct polynomial**) safe_malloc(n * n * sizeof(struct polynomial*));
        /* We don't want to actually initialize the pointers themselves, just the array containing them */
    }

    /* Allocate memory for the array B of elements on the diagonal of the final diagonal matrix */
    struct polynomial** B = (struct polynomial**) safe_malloc(A_length * sizeof(struct polynomial*));
    /* We don't want to initialize the pointers in B for the same reason */

    /* Compute s_array */
    int* s_array = (int*) safe_malloc(n * sizeof(int));
    int max;
    for (int c = 0; c < n; c++) {
        max = 0;
        for (int r = 0; r < n; r++) {
            if (MATRIX_ELEMENT(F, r, c).degree > max) max = MATRIX_ELEMENT(F, r, c);
        }
        s_array[c] = max;
    }

    /* Call the recurser function */
    polynomial_matrix_inverse_recurser(F, s_array, A, 0, B);

    /* Compute big matrix product A_1 A_2 ... A_last B^-1 = F^-1 */
    /* code later */

    /* Free memory for A */
    for (int i = 0; i < A_length; i++) {
        for (int i = 0; i < n * n; i++) {
            safe_free(A[i].data[i]->coeffs);
            safe_free(A[i].data[i]);
        }

        safe_free(A[i].data);
    }
    safe_free(A);

    /* Free memory for B */
    for (int i = 0; i < n; i++) {
        safe_free(B[i]->coeffs);
        safe_free(B[i]);
    }
    safe_free(B);

    /* s_array is already freed by the recurser algorithm */

}
