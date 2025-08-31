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

/* Function to multiply two bivariate polynomials */
struct bivariate_polynomial multiply_bivariate_polynomials(struct bivariate_polynomial P, struct bivariate_polynomial Q) {
	struct bivariate_polynomial product;
	product.lowest_degree_1 = P.lowest_degree_1 + Q.lowest_degree_1;
	product.lowest_degree_2 = P.lowest_degree_2 + Q.lowest_degree_2;
	product.highest_degree_1 = P.highest_degree_1 + Q.highest_degree_1;
	product.highest_degree_2 = P.highest_degree_2 + Q.highest_degree_2;
	product.coeffs = make_int_matrix(MAX_BIVARIATE_POLYNOMIAL_SIZE, MAX_BIVARIATE_POLYNOMIAL_SIZE);
	for (int i = product.lowest_degree_1; i <= product.highest_degree_1; i++) {
		for (int j = product.lowest_degree_2; j <= product.highest_degree_2; j++) {
			MATRIX_ELEMENT(product.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT) = 0;
		}
	}
    /*Every term of P is multiplied by every term of Q*/
	for (int i = P.lowest_degree_1; i <= P.highest_degree_1; i++) {
		for (int j = P.lowest_degree_2; j <= P.highest_degree_2; j++) {
			for (int k = Q.lowest_degree_1; k <= Q.highest_degree_1; k++) {
				for (int l = Q.lowest_degree_2; l <= Q.highest_degree_2; l++) {
					MATRIX_ELEMENT(product.coeffs, i + k + DEGREE_SHIFT, j + l + DEGREE_SHIFT) +=
						MATRIX_ELEMENT(P.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT) *
						MATRIX_ELEMENT(Q.coeffs, k + DEGREE_SHIFT, l + DEGREE_SHIFT);
				}
			}
		}
	}
	return product;
}
