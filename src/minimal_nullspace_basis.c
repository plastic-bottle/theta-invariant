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

void minimal_nullspace_basis(struct polynomial_matrix* F, struct polynomial_matrix* nullspace_basis)
{
    int m = F->rows;
    int n = F->cols;
    
    /* We will instead use an upper bound for s by finding the largest element */
    int s = 0;
    for (int r = 0; r < F->rows; r++) {
        for (int c = 0; c < F->cols; c++) {
            if (MATRIX_ELEMENT(F, r, c).degree > s) s = MATRIX_ELEMENT(F, r, c).degree;
        }
    }

    struct polynomial_matrix* P = (void*) 0;

    right_order_basis(F, ORDER_CONST*s, P);

    /* Product of F and P */
    struct polynomial_matrix* FP = (void*) 0;
    FP = mnb_fast_multiplication(F, P);

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

        nullspace_basis = P1;

        delete_polynomial_matrix(P);
        delete_polynomial_matrix(FP);

        return;
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
    minimal_nullspace_basis(G1, N1);

    struct polynomial_matrix* N2 = (void*) 0;
    minimal_nullspace_basis(mnb_fast_multiplication(G2, N1), N2);

    struct polynomial_matrix* Q = mnb_fast_multiplication(N1, N2);

    /* Recall that P2 has n rows */
    struct polynomial_matrix* P2Q = make_polynomial_matrix(n, Q->cols);
    /* Write P2Q calculation using the method from page 370 -----------------------------------------------------------------------*/
    


    /* Note that this also handles the case where F and nullspace_basis point to the same thing */
    if (nullspace_basis != NULL) {
        delete_polynomial_matrix(nullspace_basis);
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

    delete_polynomial_matrix(P);
    delete_polynomial_matrix(FP);
    delete_polynomial_matrix(G1);
    delete_polynomial_matrix(G2);
    delete_polynomial_matrix(N1);
    delete_polynomial_matrix(N2);
    delete_polynomial_matrix(P2Q);

    return;
}
