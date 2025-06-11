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

/*Function which uses synthetic division to return the quotient when a polynomial P is divided by
x - a*/
struct double_polynomial synthetic_division(struct double_polynomial P, double a) {
	struct double_polynomial quotient;
	quotient.degree = P.degree - 1;
	quotient.coeffs = (double*)safe_malloc(quotient.degree * sizeof(double));
	double next_coefficient = P.coeffs[P.degree];
	for (int index = P.degree - 1; index >= 0; index--) {
		quotient.coeffs[index] = next_coefficient;
		next_coefficient = next_coefficient * a + P.coeffs[index];
	}
	return quotient;
}
