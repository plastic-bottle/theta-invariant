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

 /*Function which returns the quotient when a bivariate polynomial is divided by T2 - 1, as needed for the theta algorithm*/
struct bivariate_polynomial theta_synthetic_division(struct bivariate_polynomial P) {
	struct bivariate_polynomial quotient;
	quotient.lowest_degree_1 = P.lowest_degree_1;
	quotient.lowest_degree_2 = P.lowest_degree_2;
	quotient.highest_degree_1 = P.highest_degree_1;
	quotient.highest_degree_2 = P.highest_degree_2 - 1;

	quotient.coeffs = make_int_matrix(MAX_BIVARIATE_POLYNOMIAL_SIZE, MAX_BIVARIATE_POLYNOMIAL_SIZE);
	THETA_INT* current_coefficient = (int*)safe_calloc((P.highest_degree_1 - P.lowest_degree_1 + 1) * sizeof(THETA_INT));
	for (int i = P.highest_degree_2; i > P.lowest_degree_2; i--) {
		for (int j = P.lowest_degree_1; j <= P.highest_degree_1; j++) {
			current_coefficient[j - P.lowest_degree_1] += MATRIX_ELEMENT(P.coeffs, j + DEGREE_SHIFT, i + DEGREE_SHIFT);
			MATRIX_ELEMENT(quotient.coeffs, j + DEGREE_SHIFT, i - 1 + DEGREE_SHIFT) = current_coefficient[j - P.lowest_degree_1];
		}
	}
	return quotient;
}