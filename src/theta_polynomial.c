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

	/*First, note down the sign of every crossing and compute the writhe of K*/
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

	/*Fill matrix A with its entries, scaled by a factor of T*/
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
	
	/*Find determinant and adjugate of A*/
	struct polynomial A_determinant = polynomial_matrix_determinant(A);
	struct polynomial_matrix* A_adjugate = polynomial_matrix_adjugate(A, A_determinant);

	int alexander_degree_shift = - 2 * n - 1 - (writhe + total_rotation) / 2;
	int adjugate_degree_shift = -2 * n - (writhe + total_rotation) / 2;

	struct bivariate_polynomial theta = initialize_bivariate_polynomial();
	int delta_123_coefficient = (writhe - total_rotation) / 2;
	struct bivariate_polynomial delta12_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial delta1_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial delta2_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial delta3_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial T2_minus_1_numerator = initialize_bivariate_polynomial();
	struct bivariate_polynomial T2_minus_1_delta12_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial T2_minus_1_delta1_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial T2_minus_1_delta2_coefficient = initialize_bivariate_polynomial();
	struct bivariate_polynomial T2_minus_1_delta3_coefficient = initialize_bivariate_polynomial();

	int one_coeffs[] = {1};
	struct polynomial one_polynomial = make_polynomial(0, one_coeffs);
	struct bivariate_polynomial one_bivariate = copy_polynomial(one_polynomial, 1, 0);
	struct bivariate_polynomial T1 = copy_polynomial(one_polynomial, 1, 1);
	struct bivariate_polynomial T1_inverse = copy_polynomial(one_polynomial, 1, -1);
	struct bivariate_polynomial T2 = copy_polynomial(one_polynomial, 2, 1);
	struct bivariate_polynomial T2_inverse = copy_polynomial(one_polynomial, 2, -1);
	struct bivariate_polynomial T3 = copy_polynomial(one_polynomial, 3, 1);
	struct bivariate_polynomial T3_inverse = copy_polynomial(one_polynomial, 3, -1);

	for (int crossing = 0; crossing < n; crossing++) {
		int sign = crossing_signs[crossing];
		int incoming_over;
		int incoming_under = K->crossings[crossing].data[0];
		if (sign == 1) {
			incoming_over = K->crossings[crossing].data[1];
		}
		else {
			incoming_over = K->crossings[crossing].data[3];
		}
		struct bivariate_polynomial current_term;
		struct bivariate_polynomial factor;
		struct bivariate_polynomial common_factor;

		current_term =  copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 3, adjugate_degree_shift);
		current_term = scale_bivariate_polynomial(current_term, -sign);
		delta12_coefficient = add_bivariate_polynomials(delta12_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 1, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, sign);
		delta3_coefficient = add_bivariate_polynomials(delta3_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 3, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, -sign);
		delta1_coefficient = add_bivariate_polynomials(delta1_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 3, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		factor = add_bivariate_polynomials(factor, scale_bivariate_polynomial(one_bivariate, -1));
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, -sign);
		delta1_coefficient = add_bivariate_polynomials(delta1_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 2, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 3, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		if (sign == 1)
			factor = T3;
		else
			factor = T3_inverse;
		factor = add_bivariate_polynomials(factor, scale_bivariate_polynomial(one_bivariate, -1));
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, sign);
		delta1_coefficient = add_bivariate_polynomials(delta1_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 1, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, -sign);
		delta3_coefficient = add_bivariate_polynomials(delta3_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 3, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, 2 * sign);
		delta1_coefficient = add_bivariate_polynomials(delta1_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 1, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 3, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, sign);
		delta2_coefficient = add_bivariate_polynomials(delta2_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 2, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 3, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, -sign);
		delta1_coefficient = add_bivariate_polynomials(delta1_coefficient, current_term);

		if (sign == 1)
			common_factor = T1;
		else
			common_factor = T1_inverse;
		common_factor = add_bivariate_polynomials(common_factor, scale_bivariate_polynomial(one_bivariate, -1));
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		common_factor = multiply_bivariate_polynomials(common_factor, factor);
		common_factor = scale_bivariate_polynomial(common_factor, sign);
		if (sign == -1) {
			common_factor = multiply_bivariate_polynomials(common_factor, T2);
			common_factor = scale_bivariate_polynomial(common_factor, -1);
		}

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 3, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 1, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		T2_minus_1_delta2_coefficient = add_bivariate_polynomials(T2_minus_1_delta2_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 2, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 1, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		current_term = scale_bivariate_polynomial(current_term, -1);
		T2_minus_1_delta3_coefficient = add_bivariate_polynomials(T2_minus_1_delta3_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 1, adjugate_degree_shift);
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		current_term = multiply_bivariate_polynomials(current_term, factor);
		T2_minus_1_delta3_coefficient = add_bivariate_polynomials(T2_minus_1_delta3_coefficient, current_term);

		if (sign == 1)
			common_factor = T3;
		else
			common_factor = T3_inverse;
		common_factor = add_bivariate_polynomials(common_factor, scale_bivariate_polynomial(one_bivariate, -1));
		factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 3, adjugate_degree_shift);
		common_factor = multiply_bivariate_polynomials(common_factor, factor);
		common_factor = scale_bivariate_polynomial(common_factor, sign);
		if (sign == -1) {
			common_factor = multiply_bivariate_polynomials(common_factor, T2);
			common_factor = scale_bivariate_polynomial(common_factor, -1);
		}

		current_term = one_bivariate;
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		T2_minus_1_delta12_coefficient = add_bivariate_polynomials(T2_minus_1_delta12_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_over - 1), 1, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, -1);
		T2_minus_1_delta2_coefficient = add_bivariate_polynomials(T2_minus_1_delta2_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_over - 1, incoming_under - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		T2_minus_1_delta1_coefficient = add_bivariate_polynomials(T2_minus_1_delta1_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_under - 1), 2, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		factor = add_bivariate_polynomials(factor, scale_bivariate_polynomial(one_bivariate, -2));
		current_term = multiply_bivariate_polynomials(current_term, factor);
		T2_minus_1_delta1_coefficient = add_bivariate_polynomials(T2_minus_1_delta1_coefficient, current_term);

		current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, incoming_over - 1), 1, adjugate_degree_shift);
		current_term = multiply_bivariate_polynomials(current_term, common_factor);
		if (sign == 1)
			factor = T1;
		else
			factor = T1_inverse;
		factor = add_bivariate_polynomials(factor, scale_bivariate_polynomial(one_bivariate, -1));
		current_term = multiply_bivariate_polynomials(current_term, factor);
		if (sign == 1)
			factor = T2;
		else
			factor = T2_inverse;
		factor = add_bivariate_polynomials(factor, one_bivariate);
		current_term = multiply_bivariate_polynomials(current_term, factor);
		current_term = scale_bivariate_polynomial(current_term, -1);
		T2_minus_1_delta2_coefficient = add_bivariate_polynomials(T2_minus_1_delta2_coefficient, current_term);

		for (int other_crossing = 0; other_crossing < n; other_crossing++) {
			int other_sign = crossing_signs[other_crossing];
			int other_incoming_over;
			int other_incoming_under = K->crossings[other_crossing].data[0];
			if (other_sign == 1) {
				other_incoming_over = K->crossings[other_crossing].data[1];
			}
			else {
				other_incoming_over = K->crossings[other_crossing].data[3];
			}

			common_factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, other_incoming_under - 1, incoming_over - 1), 1, adjugate_degree_shift);
			factor = copy_polynomial(MATRIX_ELEMENT(A_adjugate, incoming_under - 1, other_incoming_over - 1), 3, adjugate_degree_shift);
			common_factor = multiply_bivariate_polynomials(common_factor, factor);
			if (sign == 1)
				factor = T1;
			else
				factor = T1_inverse;
			factor = add_bivariate_polynomials(factor, scale_bivariate_polynomial(one_bivariate, -1));
			common_factor = multiply_bivariate_polynomials(common_factor, factor);
			if (other_sign == 1)
				factor = T3;
			else
				factor = T3_inverse;
			factor = add_bivariate_polynomials(factor, scale_bivariate_polynomial(one_bivariate, -1));
			common_factor = multiply_bivariate_polynomials(common_factor, factor);
			common_factor = scale_bivariate_polynomial(common_factor, other_sign);
			if (other_sign == -1) {
				common_factor = multiply_bivariate_polynomials(common_factor, T2);
				common_factor = scale_bivariate_polynomial(common_factor, -1);
			}

			current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, other_incoming_over - 1, incoming_over - 1), 2, adjugate_degree_shift);
			current_term = multiply_bivariate_polynomials(current_term, common_factor);
			if (sign == 1)
				factor = T2;
			else
				factor = T2_inverse;
			current_term = multiply_bivariate_polynomials(current_term, factor);
			T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, current_term);

			current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, other_incoming_under - 1, incoming_under - 1), 2, adjugate_degree_shift);
			current_term = multiply_bivariate_polynomials(current_term, common_factor);
			T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, current_term);

			current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, other_incoming_under - 1, incoming_over - 1), 2, adjugate_degree_shift);
			current_term = multiply_bivariate_polynomials(current_term, common_factor);
			if (sign == 1)
				factor = T2;
			else
				factor = T2_inverse;
			current_term = multiply_bivariate_polynomials(current_term, factor);
			current_term = scale_bivariate_polynomial(current_term, -1);
			T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, current_term);

			current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, other_incoming_over - 1, incoming_under - 1), 2, adjugate_degree_shift);
			current_term = multiply_bivariate_polynomials(current_term, common_factor);
			current_term = scale_bivariate_polynomial(current_term, -1);
			T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, current_term);
		}
	}

	for (int strand = 1; strand <= 2 * n + 1; strand++) {
		struct bivariate_polynomial current_term = copy_polynomial(MATRIX_ELEMENT(A_adjugate, strand - 1, strand - 1), 3, adjugate_degree_shift);
		current_term = scale_bivariate_polynomial(current_term, rotation[strand]);
		delta12_coefficient = add_bivariate_polynomials(delta12_coefficient, current_term);
	}

	struct bivariate_polynomial delta1 = copy_polynomial(A_determinant, 1, alexander_degree_shift);
	struct bivariate_polynomial delta2 = copy_polynomial(A_determinant, 2, alexander_degree_shift);
	struct bivariate_polynomial delta3 = copy_polynomial(A_determinant, 3, alexander_degree_shift);
	struct bivariate_polynomial delta12 = multiply_bivariate_polynomials(delta1, delta2);	
	struct bivariate_polynomial delta123 = multiply_bivariate_polynomials(delta12, delta3);	

	struct bivariate_polynomial delta1_product = multiply_bivariate_polynomials(delta1, delta1_coefficient);
	theta = add_bivariate_polynomials(theta, delta1_product);
	struct bivariate_polynomial delta2_product = multiply_bivariate_polynomials(delta2, delta2_coefficient);
	theta = add_bivariate_polynomials(theta, delta2_product);
	struct bivariate_polynomial delta3_product = multiply_bivariate_polynomials(delta3, delta3_coefficient);
	theta = add_bivariate_polynomials(theta, delta3_product);
	struct bivariate_polynomial delta12_product = multiply_bivariate_polynomials(delta12, delta12_coefficient);
	theta = add_bivariate_polynomials(theta, delta12_product);
	struct bivariate_polynomial delta123_product = scale_bivariate_polynomial(delta123, delta_123_coefficient);
	theta = add_bivariate_polynomials(theta, delta123_product);

	struct bivariate_polynomial T2_minus_1_delta1_product = multiply_bivariate_polynomials(delta1, T2_minus_1_delta1_coefficient);
	T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, T2_minus_1_delta1_product);
	struct bivariate_polynomial T2_minus_1_delta2_product = multiply_bivariate_polynomials(delta2, T2_minus_1_delta2_coefficient);
	T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, T2_minus_1_delta2_product);
	struct bivariate_polynomial T2_minus_1_delta3_product = multiply_bivariate_polynomials(delta3, T2_minus_1_delta3_coefficient);
	T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, T2_minus_1_delta3_product);
	struct bivariate_polynomial T2_minus_1_delta12_product = multiply_bivariate_polynomials(delta12, T2_minus_1_delta12_coefficient);
	T2_minus_1_numerator = add_bivariate_polynomials(T2_minus_1_numerator, T2_minus_1_delta12_product);
	struct bivariate_polynomial T2_minus_1_quotient = theta_synthetic_division(T2_minus_1_numerator);
	theta = add_bivariate_polynomials(theta, T2_minus_1_quotient);

	return theta;
}