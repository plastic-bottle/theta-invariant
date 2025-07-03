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

/* Helper function for polynomial_matrix_inverse */
/* DO NOT CALL THIS OUTSIDE OF polynomial_matrix_inverse */
/* F is the input matrix to invert, A is the pointer to the particular A_k that this function call is concerned with, B is the resulting diagonal matrix */
/* This function will handle freeing the memory in F if needed */
void polynomial_matrix_inverse_recurser(struct polynomial_split_matrix* const F, int* s_array, struct polynomial_split_matrix* const A, struct polynomial** const B)
{
    /* */
    /* do we need to check if F is zero or something dumb? */
    /* */

    if (F->rows == 1) {
        /* The element in B will point to the same memory as the first row pointer in F->data (which points to the first element of the first row ) */
        /* Don't free F since we will use that memory in B */
        B[0] = F->data[0];

        return;
    }

    

}
