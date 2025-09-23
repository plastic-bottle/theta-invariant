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



/* Puts results in order_basis and s_col_degs */
/* If the input and result arguments are identical, then old memory will be freed and new memory will be allocated */
void right_order_basis(struct polynomial_matrix* F, int order, int* s_array, struct polynomial_matrix* order_basis, int* s_col_degs)
{

}

struct polynomial_matrix* shift_and_truncate(struct polynomial_matrix* A, int shift, int maxDeg)
{
    struct polynomial_matrix* B = make_polynomial_matrix(A->rows, A->cols);
    for (int i=0; i<A->rows; i++) {
        for (int j=0; j<A->cols; j++) {
            int newDeg = MATRIX_ELEMENT(A, i, j).degree + shift;
            if (newDeg<0) newDeg=0;
            if (newDeg>maxDeg) newDeg=maxDeg;

            struct polynomial temp;
            temp.degree = newDeg;
            temp.coeffs = (int*) safe_malloc((temp.degree + 1) * sizeof(int));
            MATRIX_ELEMENT(B, i, j) = temp;

            for (int k=0; k<=MATRIX_ELEMENT(A, i, j).degree; k++) {
                int idx = k + shift;
                if(idx>=0 && idx<=maxDeg) MATRIX_ELEMENT(B, i, j).coeffs[idx] = MATRIX_ELEMENT(A, i, j).coeffs[k];
            }
        }
    }
    return B;
}