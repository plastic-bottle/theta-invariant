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

 /* Function to find the determinant of a polynomial matrix using the Bareiss algorithm */
struct polynomial polynomial_matrix_determinant(struct polynomial_matrix* A) {
	size_t n = A->rows;
	struct polynomial_matrix* A_copy = make_polynomial_matrix(n, n);
	A_copy->rows = A_copy->cols = n;
	A_copy->data = (struct polynomial*)safe_malloc(n * n * sizeof(struct polynomial));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			MATRIX_ELEMENT(A_copy, i, j) = MATRIX_ELEMENT(A, i, j);
		}
	}
	int one_coeffs[] = { 1 };
	struct polynomial one_polynomial = make_polynomial(0, one_coeffs);
	struct polynomial denominator = one_polynomial;

	for (size_t k = 0; k < n - 1; k++) {
		for (size_t i = k + 1; i < n; i++) {
			for (size_t j = k + 1; j < n; j++) {
				struct polynomial first_product = multiply_polynomials(MATRIX_ELEMENT(A_copy, i, j), MATRIX_ELEMENT(A_copy, k, k));
				struct polynomial second_product = multiply_polynomials(MATRIX_ELEMENT(A_copy, i, k), MATRIX_ELEMENT(A_copy, k, j));
				struct polynomial difference = subtract_polynomials(first_product, second_product);
				MATRIX_ELEMENT(A_copy, i, j) = divide_polynomials(difference, denominator);
			}
			MATRIX_ELEMENT(A_copy, i, k) = initialize_polynomial();
		}
		denominator = MATRIX_ELEMENT(A_copy, k, k);
	}
	return MATRIX_ELEMENT(A_copy, n - 1, n - 1);
}