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

/* Base case for right order basis computation */
struct polynomial_matrix* right_order_basis_base_case(struct int_matrix* A)
{
    int rows = A->rows, cols = A->cols;
    struct int_matrix* B = make_int_matrix(cols, cols);

    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++) {
            MATRIX_ELEMENT(B, i, j) = MATRIX_ELEMENT(A, i, j);
        }
    }

    int* one = {1};

    struct polynomial_matrix* N = make_polynomial_matrix(cols, cols);
    for(int i=0; i<cols; i++){
        MATRIX_ELEMENT(N, i, i) = make_polynomial(0, one);
    }

    int *pivotCols = (int*)calloc(cols, sizeof(int));

    for (int c=0; c<cols; c++) {
        int pivotRow=-1;
        for (int r=c; r<cols; r++) if (MATRIX_ELEMENT(B, r, c) != 0) {pivotRow=r; break;}
        if (pivotRow==-1) continue;
        pivotCols[c]=1;
        if (pivotRow!=c) for (int j=0;j<cols;j++) SWAP(MATRIX_ELEMENT(B, pivotRow, j),MATRIX_ELEMENT(B, c, j), int);

        for (int j=0; j<cols; j++) {
            if (j == c) continue;
            int f1 = MATRIX_ELEMENT(B, c, c);
            int f2 = MATRIX_ELEMENT(B, c, j);
            for (int i=0; i<cols; i++){
                MATRIX_ELEMENT(B, i, j) = f1 * MATRIX_ELEMENT(B, i, j) - f2 * MATRIX_ELEMENT(B, i, c);
                MATRIX_ELEMENT(N, i, j).coeffs[0] = f1 * MATRIX_ELEMENT(N, i, j).coeffs[0] - f2 * MATRIX_ELEMENT(N, i, c).coeffs[0];
            }
        }
    }

    struct polynomial p, old;
    for (int j=0; j<cols; j++){
        for (int i=0; i<cols; i++){
            old = MATRIX_ELEMENT(N, i, j);
            if (pivotCols[j]) {
                old.degree++;
                for (int k=0; k<=old.degree; k++) old.coeffs[k+1] = old.coeffs[k];
                old.coeffs[0] = 0;
            }
        }
    }

    freeMatrix(B);
    free(pivotCols);
    return N;
}