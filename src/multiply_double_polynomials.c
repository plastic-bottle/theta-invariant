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

/*Function which returns the product of two double polynomials*/
struct double_polynomial multiply_double_polynomials(struct double_polynomial P, struct double_polynomial Q) {
	struct double_polynomial product;
	product.degree = P.degree + Q.degree;
	product.coeffs = (int*)safe_malloc(MAX_POLYNOMIAL_SIZE * sizeof(double));
	for (int index = 0; index <= product.degree; index++) 
		product.coeffs[index] = 0;
    /*Each term from P is multiplied with every other term from Q*/
	for (int P_index = 0; P_index <= P.degree; P_index++) 
		for (int Q_index = 0; Q_index <= Q.degree; Q_index++) 
			product.coeffs[P_index + Q_index] += P.coeffs[P_index] * Q.coeffs[Q_index];
	return product;
}