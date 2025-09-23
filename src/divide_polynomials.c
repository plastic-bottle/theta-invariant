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

/* Calculates the quotient P/Q, throws away remainder */
struct polynomial divide_polynomials(struct polynomial P, struct polynomial Q)
{
    if (P.degree == 0 && P.coeffs[0] == 0) {
        return initialize_polynomial();
    }
    
    struct polynomial quotient;
	quotient.degree = P.degree - Q.degree;
	quotient.coeffs = (int*) safe_malloc((quotient.degree + 1) * sizeof(int));

    int current_coeff;
	for (int result_index = (int) quotient.degree; result_index >= 0; result_index--) {
        current_coeff = P.coeffs[result_index + (int) Q.degree] / Q.coeffs[(int) Q.degree];
        quotient.coeffs[result_index] = current_coeff;
        for (int intermediate_index = (int) Q.degree; intermediate_index >= 0; intermediate_index--) {
		    P.coeffs[result_index + intermediate_index] -= current_coeff * Q.coeffs[intermediate_index];
        }
	}

	return quotient;
}