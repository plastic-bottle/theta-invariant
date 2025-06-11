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

/*Function which, when given the maximum degree n for a Lagrange interpolating polynomials and n + 1
inputs for x, returns the product of all of the x - x_i*/
struct double_polynomial lagrange_product(int max_degree, double* inputs) {
	struct double_polynomial product;
	product.degree = 0;
	product.coeffs = (double*)safe_malloc(sizeof(double));
	product.coeffs[0] = 1;
	for (int index = 0; index <= max_degree; index++) {
		double factor_coefficients[] = { -inputs[index], 1};
		struct double_polynomial factor = make_double_polynomial(1, factor_coefficients);
		product = multiply_double_polynomials(product, factor);
	}
	return product;
}
