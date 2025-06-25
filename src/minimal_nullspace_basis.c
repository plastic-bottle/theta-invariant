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

/* Requires F to be s-reduced */
/* Puts results in nullspace_basis and s_col_degs */
/* If any input and result arguments are identical, then old memory will be freed and new memory will be allocated */
void minimal_nullspace_basis(struct polynomial_matrix* F, int* s_array, struct polynomial_matrix* nullspace_basis, int* s_col_degs)
{
    int m = F->rows;
    int n = F->cols;
    
    /* We will instead use an upper bound for s by finding the largest element */
    int s = 0;
    for (int i = 0; i < n; i++) {
        if (s_array[i] > s) s = s_array[i];
    }

    struct polynomial_matrix* P = make_polynomial_matrix(n, n);
    /* s-column degrees of P */
    int* b_array = (int*) safe_malloc(n * sizeof(int));

    right_order_basis(F, 3*s, s_array, P, b_array);

    /* Product of F and P */
    struct polynomial_matrix* FP = make_polynomial_matrix(m, n);
    FP = mnb_fast_multiplication(F, P);

    /* Assuming that all of P1 will come at the left of P - should verify this assumption later */
    int P1_cols = 0;
    bool end_loop = false;
    /* The number of cols in P1 is bounded by the number of cols in the nullspace basis, which is in turn bounded by the number of rows in F */
    for (int c = 0; c < m && end_loop == false; c++) {
        for (int r = 0; r < m && end_loop == false; r++) {
            if (MATRIX_ELEMENT(FP, r, c).degree != 0 || MATRIX_ELEMENT(FP, r, c).coeffs[0] != 0) end_loop = true;
        }
        
        if (end_loop == false) {
            P1_cols++;
        }
    }

    /* Copy everything in P1 over to its own matrix */
    /* This can be optimized out if need be, but you would need to be able to selectively free parts of P */
    struct polynomial_matrix* P1 = make_polynomial_matrix(n, P1_cols);
    for (int c = 0; c < P1->rows; c++) {
        for (int r = 0 ; r < P1->cols; r++) {
            MATRIX_ELEMENT(P1, r, c) = make_polynomial(MATRIX_ELEMENT(P, r, c).degree, MATRIX_ELEMENT(P, r, c).coeffs);
        }
    }

    P1_s_col_degs = (int*) safe_malloc(P1_cols * sizeof(int));
    for (int i = 0; i < P1_cols; i++) {
        s_col_degs[i] = b_array[i];
    }

    if (m == 1) {
        /* Note that this also handles the case where F and nullspace_basis point to the same thing */
        if (nullspace_basis != NULL) {
            delete_polynomial_matrix(nullspace_basis);
        }

        /* Note that this also handles the case where s_array and s_col_degs point to the same thing */
        if (s_col_degs != NULL) {
            safe_free(s_col_degs);
        }

        nullspace_basis = P1;

        s_col_degs = P1_s_col_degs;

        delete_polynomial_matrix(P);
        delete_polynomial_matrix(FP);
        safe_free(b_array);

        return;
    }

    /* t_array contains the downshifted s-column degrees of P2 */
    int* t_array = (int*) safe_malloc((n - P1_cols) * sizeof(int));
    for (int i = 0; i < (n - P1_cols); i++) {
        t_array[i] = b_array[i + P1_cols] - 3*s;
    }

    /* Copy the portion of FP representing FP2 into an array called G while dividing by x^(3s) */
    struct polynomial_matrix* G = make_polynomial_matrix(n, n - P1_cols);
    for (int c = 0; c < G->rows; c++) {
        for (int r = 0 ; r < G->cols; r++) {
            
        }
    }

}
