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

 /*Function which returns the sum of two bivariate polynomials*/
struct bivariate_polynomial add_bivariate_polynomials(struct bivariate_polynomial P, struct bivariate_polynomial Q) {
	struct bivariate_polynomial sum;
	sum.lowest_degree_1 = MIN(P.lowest_degree_1, Q.lowest_degree_1);
	sum.lowest_degree_2 = MIN(P.lowest_degree_2, Q.lowest_degree_2);
	sum.highest_degree_1 = MAX(P.highest_degree_1, Q.highest_degree_1);
	sum.highest_degree_2 = MAX(P.highest_degree_2, Q.highest_degree_2);
	sum.coeffs = make_int_matrix(MAX_BIVARIATE_POLYNOMIAL_SIZE, MAX_BIVARIATE_POLYNOMIAL_SIZE);

    /*Each coefficient of the sum is the sum of the two corresponding coefficients from P and Q*/
	for (int i = sum.lowest_degree_1; i <= sum.highest_degree_1; i++) {
		for (int j = sum.lowest_degree_2; j <= sum.highest_degree_2; j++) {
			MATRIX_ELEMENT(sum.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT) = 0;
			if (i >= P.lowest_degree_1 && i <= P.highest_degree_1 && j >= P.lowest_degree_2 && j <= P.highest_degree_2) {
				MATRIX_ELEMENT(sum.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT) += MATRIX_ELEMENT(P.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT);
			}
			if (i >= Q.lowest_degree_1 && i <= Q.highest_degree_1 && j >= Q.lowest_degree_2 && j <= Q.highest_degree_2) {
				MATRIX_ELEMENT(sum.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT) += MATRIX_ELEMENT(Q.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT);
			}
		}
	}
	return sum;
}
