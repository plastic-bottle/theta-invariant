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

 /*Function which computes the theta polynomial of a knot K*/
struct bivariate_polynomial theta_polynomial(const struct knot* const K) {
	int n = K->number_of_crossings;

	int* crossing_signs = (int*)safe_malloc(n * sizeof(int));
	int writhe = 0;
	for (int crossing = 0; crossing < n; crossing++)
		if ((K->crossings[crossing].data[3] - K->crossings[crossing].data[1] + 2 * n) % (2 * n) == 1) {
			writhe++;
			crossing_signs[crossing] = 1;
		} else {
			writhe--;
			crossing_signs[crossing] = -1;
		}
	int* rotation = rotation_numbers(K);
	int total_rotation = 0;
	for (int strand = 1; strand <= 2 * n + 1; strand++)
		total_rotation += rotation[strand];
	int alexander_degree_shift = (-writhe - total_rotation) / 2;

	int T_coeffs[] = { 0, 1 };
	struct polynomial T = make_polynomial(1, T_coeffs);
	struct polynomial_matrix* A = (struct polynomial_matrix*)safe_malloc(sizeof(struct polynomial_matrix));
	A->rows = A->cols = 2 * n + 1;
	A->data = (struct polynomial*)safe_malloc(A->rows * A->cols * sizeof(struct polynomial));
	for (int i = 0; i < A->rows; i++)	
		for (int j = 0; j < A->cols; j++) {
			MATRIX_ELEMENT(A, i, j) = initialize_polynomial();
			if (i == j)
				MATRIX_ELEMENT(A, i, j) = T;
		}
	int negative_one_coeffs[] = { -1 };
	struct polynomial negative_one = make_polynomial(0, negative_one_coeffs);
	int one_minus_T_coeffs[] = { 1, -1 };
	struct polynomial one_minus_T = make_polynomial(1, one_minus_T_coeffs);
	int minus_T_coeffs[] = { 0, -1 };
	struct polynomial minus_T = make_polynomial(1, minus_T_coeffs);
	int minus_T_squared_coeffs[] = { 0, 0, -1 };
	struct polynomial minus_T_squared = make_polynomial(2, minus_T_squared_coeffs);
	int T_squared_minus_T_coeffs[] = { 0, -1, 1 };
	struct polynomial T_squared_minus_T = make_polynomial(2, T_squared_minus_T_coeffs);
	for (int crossing = 0; crossing < n; crossing++) {
		if (crossing_signs[crossing] == 1) {
			int incoming_over = K->crossings[crossing].data[1];
			int incoming_under = K->crossings[crossing].data[0];
			MATRIX_ELEMENT(A, incoming_over - 1, incoming_over) = minus_T_squared;
			MATRIX_ELEMENT(A, incoming_over - 1, incoming_under) = T_squared_minus_T;
			MATRIX_ELEMENT(A, incoming_under - 1, incoming_under) = minus_T;
		} else {
			int incoming_over = K->crossings[crossing].data[3];
			int incoming_under = K->crossings[crossing].data[0];
			MATRIX_ELEMENT(A, incoming_over - 1, incoming_over) = negative_one;	
			MATRIX_ELEMENT(A, incoming_over - 1, incoming_under) = one_minus_T;
			MATRIX_ELEMENT(A, incoming_under - 1, incoming_under) = minus_T;
		}
	}
	
	//find adjugate of A and its determinant to get Alexander polynomial
	struct polynomial_matrix* A_adjugate = (struct polynomial_matrix*)safe_malloc(sizeof(struct polynomial_matrix));
	A_adjugate->rows = A->rows;
	A_adjugate->cols = A->cols;
	A_adjugate->data = (struct polynomial*)safe_malloc(A->rows * A->cols * sizeof(struct polynomial));
	struct polynomial shifted_alexander_polynomial;

    //compute theta
	struct bivariate_polynomial theta = initialize_bivariate_polynomial();

	return theta;
}