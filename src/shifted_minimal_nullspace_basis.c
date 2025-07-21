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
void shifted_minimal_nullspace_basis(struct polynomial_matrix* F, int* s_array, struct polynomial_matrix* nullspace_basis, int* s_col_degs)
{
    int m = F->rows;
    int n = F->cols;
    
    /* We will instead use an upper bound for s by finding the largest element */
    int s = 0;
    for (int i = 0; i < n; i++) {
        if (s_array[i] > s) s = s_array[i];
    }

    struct polynomial_matrix* P = (void*) 0;
    /* s-column degrees of P */
    int* b_array = (int*) safe_malloc(n * sizeof(int));

    right_order_basis(F, ORDER_CONST*s, s_array, P, b_array);

    /* Product of F and P */
    struct polynomial_matrix* FP = multiply_polynomial_matrices(F, P);

    /* Assuming that all of P1 will come at the left of P - should verify this assumption later */
    int P1_cols = 0;
    enum boolean end_loop = FALSE;
    /* The number of cols in P1 is bounded by the number of cols in the nullspace basis, which is in turn bounded by the number of rows in F */
    for (int c = 0; c < m && end_loop == FALSE; c++) {
        for (int r = 0; r < m && end_loop == FALSE; r++) {
            if (MATRIX_ELEMENT(FP, r, c).degree != 0 || MATRIX_ELEMENT(FP, r, c).coeffs[0] != 0) end_loop = TRUE;
        }
        
        if (end_loop == FALSE) {
            P1_cols++;
        }
    }

    if (m == 1) {
        /* Copy everything in P1 over to its own matrix */
        /* This can be optimized out if need be, but you would need to be able to selectively free parts of P */
        struct polynomial_matrix* P1 = make_polynomial_matrix(n, P1_cols);
        for (int c = 0; c < P1->cols; c++) {
            for (int r = 0 ; r < P1->rows; r++) {
                MATRIX_ELEMENT(P1, r, c) = make_polynomial(MATRIX_ELEMENT(P, r, c).degree, MATRIX_ELEMENT(P, r, c).coeffs);
            }
        }

        /* Note that this also handles the case where F and nullspace_basis point to the same thing */
        if (nullspace_basis != NULL) {
            delete_polynomial_matrix(nullspace_basis);
        }

        /* Note that this also handles the case where s_array and s_col_degs point to the same thing */
        if (s_col_degs != NULL) {
            safe_free(s_col_degs);
        }

        nullspace_basis = P1;

        s_col_degs = (int*) safe_malloc(P1_cols * sizeof(int));
        for (int i = 0; i < P1_cols; i++) {
            s_col_degs[i] = b_array[i];
        }

        delete_polynomial_matrix(P);
        delete_polynomial_matrix(FP);
        safe_free(b_array);

        return;
    }

    /* t_array contains the downshifted s-column degrees of P2 */
    int* t_array = (int*) safe_malloc((n - P1_cols) * sizeof(int));
    for (int i = 0; i < (n - P1_cols); i++) {
        t_array[i] = b_array[i + P1_cols] - ORDER_CONST*s;
    }

    /* Copy the portion of FP representing FP2 into G1 and G2 while dividing by x^(ORDER_CONST*s) */
    struct polynomial_matrix* G1 = make_polynomial_matrix((int) (m / 2), n - P1_cols);
    struct polynomial_matrix* G2 = make_polynomial_matrix(m - G1->rows, n - P1_cols);
    for (int c = 0; c < n - P1_cols; c++) {
        for (int r = 0; r < G1->rows; r++) {
            MATRIX_ELEMENT(G1, r, c) = initialize_polynomial();
            MATRIX_ELEMENT(G1, r, c).degree = MATRIX_ELEMENT(FP, r, c + P1_cols).degree - ORDER_CONST*s;
            for (int i = 0; i <= MATRIX_ELEMENT(G1, r, c).degree; i++) {
                MATRIX_ELEMENT(G1, r, c).coeffs[i] = MATRIX_ELEMENT(FP, r, c + P1_cols).coeffs[i + ORDER_CONST*s];
            }
        }

        for (int r = 0; r < G2->rows; r++) {
            MATRIX_ELEMENT(G2, r, c) = initialize_polynomial();
            MATRIX_ELEMENT(G2, r, c).degree = MATRIX_ELEMENT(FP, r + G1->rows, c + P1_cols).degree - ORDER_CONST*s;
            for (int i = 0; i <= MATRIX_ELEMENT(G2, r, c).degree; i++) {
                MATRIX_ELEMENT(G2, r, c).coeffs[i] = MATRIX_ELEMENT(FP, r + G1->rows, c + P1_cols).coeffs[i + ORDER_CONST*s];
            }
        }
    }

    struct polynomial_matrix* N1 = (void*) 0;
    int* u_array = (void*) 0;
    shifted_minimal_nullspace_basis(G1, t_array, N1, u_array);

    struct polynomial_matrix* G2N1 = multiply_polynomial_matrices(G2, N1);

    struct polynomial_matrix* N2 = (void*) 0;
    int* v_array = (void*) 0;
    shifted_minimal_nullspace_basis(G2N1, u_array, N2, v_array);

    struct polynomial_matrix* Q = multiply_polynomial_matrices(N1, N2);

    struct polynomial_matrix* P2Q = multiply_polynomial_matrices(P2, Q);

    /* Note that this also handles the case where F and nullspace_basis point to the same thing */
    if (nullspace_basis != NULL) {
        delete_polynomial_matrix(nullspace_basis);
    }

    /* Note that this also handles the case where s_array and s_col_degs point to the same thing */
    if (s_col_degs != NULL) {
        safe_free(s_col_degs);
    }

    nullspace_basis = make_polynomial_matrix(n, P1_cols + P2Q->cols);
    for (int r = 0; r < nullspace_basis->rows; r++) {
        /* Copy elements from P1 over */
        for (int c = 0; c < P1_cols; c++) {
            MATRIX_ELEMENT(nullspace_basis, r, c) = make_polynomial(MATRIX_ELEMENT(P, r, c).degree, MATRIX_ELEMENT(P, r, c).coeffs);
        }
        
        /* Copy elements from P2Q over */
        for (int c = 0; c < P2Q->cols; c++) {
            MATRIX_ELEMENT(nullspace_basis, r, c + P1_cols) = make_polynomial(MATRIX_ELEMENT(P2Q, r, c).degree, MATRIX_ELEMENT(P2Q, r, c).coeffs);
        }
    }

    s_col_degs = (int*) safe_malloc((P1_cols + P2Q->cols) * sizeof(int));
    for (int i = 0; i < P1_cols; i++) {
        s_col_degs[i] = b_array[i];
    }
    for (int i = 0; i < P2Q->cols; i++) {
        /* I think that we need to shift v back by ORDER_CONST*c */
        s_col_degs[i + P1_cols] = v_array[i] + ORDER_CONST*c;
    }

    delete_polynomial_matrix(P);
    delete_polynomial_matrix(FP);
    delete_polynomial_matrix(G1);
    delete_polynomial_matrix(G2);
    delete_polynomial_matrix(N1);
    delete_polynomial_matrix(N2);
    delete_polynomial_matrix(G2N1);
    delete_polynomial_matrix(P2Q);
    safe_free(b_array);
    safe_free(t_array);
    safe_free(u_array);
    safe_free(v_array);

    return;

}
