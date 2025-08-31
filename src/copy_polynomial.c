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

 /* Function to copy a one-variable polynomial into a bivariate polynomial, depending on the variable it is to be copied to */
struct bivariate_polynomial copy_polynomial(struct polynomial P, int variable_index, int shift) {
	struct bivariate_polynomial copy;
	copy.coeffs = make_int_matrix(MAX_BIVARIATE_POLYNOMIAL_SIZE, MAX_BIVARIATE_POLYNOMIAL_SIZE);
    /*If variable index is 1, P is copied into a polynomial in T1*/
	if (variable_index == 1) {
		copy.lowest_degree_1 = shift;
		copy.highest_degree_1 = P.degree + shift;
		copy.lowest_degree_2 = copy.highest_degree_2 = 0;
		for (int i = 0; i <= P.degree; i++) {
			MATRIX_ELEMENT(copy.coeffs, i + shift + DEGREE_SHIFT, DEGREE_SHIFT) = P.coeffs[i];
		}
    /*If variable index is 2, P is copied into a polynomial in T2*/
	} else if (variable_index == 2) {
		copy.lowest_degree_1 = copy.highest_degree_1 = 0;
		copy.lowest_degree_2 = shift;
		copy.highest_degree_2 = P.degree + shift;
		for (int i = 0; i <= P.degree; i++) {
			MATRIX_ELEMENT(copy.coeffs, DEGREE_SHIFT, i + shift + DEGREE_SHIFT) = P.coeffs[i];
		}
    /*If variable index is 3, P is copied into a polynomial in T3 = T1 * T2*/
	} else if (variable_index == 3) {
		copy.lowest_degree_1 = copy.lowest_degree_2 = shift;
		copy.highest_degree_1 = copy.highest_degree_2 = P.degree + shift;
		for (int i = 0; i <= P.degree; i++) {
			for (int j = 0; j <= P.degree; j++) {
				if (i == j) {
					MATRIX_ELEMENT(copy.coeffs, i + shift + DEGREE_SHIFT, j + shift + DEGREE_SHIFT) = P.coeffs[i];
				}
				else {
					MATRIX_ELEMENT(copy.coeffs, i + shift + DEGREE_SHIFT, j + shift + DEGREE_SHIFT) = 0;
				}
			}
		}
	}

	return copy;
}
