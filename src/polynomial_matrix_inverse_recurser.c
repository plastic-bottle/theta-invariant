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
void polynomial_matrix_inverse_recurser(struct polynomial_matrix* const F, int* s_array, struct polynomial_pointer_matrix* const A, int a_start, struct polynomial** const B)
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

    int chop = F->rows / 2;

    /* Initialize F_u and F_d */
    struct polynomial_matrix* F_u = (struct polynomial_matrix*) safe_malloc(sizeof(struct polynomial_matrix));
    result->rows = chop;
    result->cols = F->cols;
    struct polynomial_matrix* F_d = (struct polynomial_matrix*) safe_malloc(sizeof(struct polynomial_matrix));
    result->rows = F->rows - chop;
    result->cols = F->cols;
    /* Transfer memory from F to F_u and F_d */
    F_u->data = F->data;
    F_d->data = F->data + chop * F->cols;

    /* Calculate the nullspace bases */
    struct polynomial_matrix* N_r = (void*) 0;
    int* N_r_cdegs = (void*) 0;
    shifted_minimal_nullspace_basis(F_u, s_array, N_r, N_r_cdegs);

    struct polynomial_matrix* N_l = (void*) 0;
    int* N_l_cdegs = (void*) 0;
    shifted_minimal_nullspace_basis(F_d, s_array, N_l, N_l_cdegs);

    /* Calculate the new entries of the working matrix after multiplication with the nullspace matrix */
    struct polynomial_matrix* R_u = make_polynomial_pointer_matrix(F_u->rows, N_l->cols);
    /* R_u = F_u * N_l */

    struct polynomial_matrix* R_d = make_polynomial_pointer_matrix(F_d->rows, N_r->cols);
    /* R_d = F_d * N_r */

    /* Update A with our nullspace matrix */
    for (int r = 0; r < F->rows; r++) {
        /* Copy elements from N_l over */
        for (int c = 0; c < F->rows-chop; c++) {
            MATRIX_ELEMENT(A, r + a_start, c + a_start) = &MATRIX_ELEMENT(N_l, r, c);
        }
        /* Copy elements from N_r over */
        for (int c = 0; c < chop; c++) {
            MATRIX_ELEMENT(A, r + a_start, c + a_start + F->rows-chop) = &MATRIX_ELEMENT(N_r, r, c);
        }
    }

    /* Free F and stuff attached to it */
    delete_polynomial_matrix(F);
    safe_free(F_u->data);
    safe_free(F_d->data);
    safe_free(F_u);
    safe_free(F_d);

    /* Free s_array */
    safe_free(s_array);

    /* Free pointers to N_r and N_l but not their contents */
    safe_free(N_r);
    safe_free(N_l);

    polynomial_matrix_inverse_recurser(R_u, N_l_cdegs, A+1, a_start, B);
    polynomial_matrix_inverse_recurser(R_d, N_r_cdegs, A+1, a_start+chop, B+chop);
    
    /* No need to free the N_l and N_r cdegs because the function call itself frees it */

}
