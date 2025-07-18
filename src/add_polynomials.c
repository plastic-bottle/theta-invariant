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

struct polynomial add_polynomials(struct polynomial P, struct polynomial Q)
{
    struct polynomial sum;
    sum.degree = MAX(P.degree, Q.degree);
    if (P.degree == Q.degree) { /* Adjust degree if leading terms cancel out */
        while (P.coeffs[sum.degree] == -1 * Q.coeffs[sum.degree] && sum.degree > 0) {
            sum.degree -= 1
        }
    }
	sum.coeffs = (int*) safe_malloc(MAX_POLYNOMIAL_SIZE * sizeof(int));

    /*Each of the coefficients of the sum is set to the sum of the corresponding coefficients from P and Q*/
	for (int index = 0; index <= sum.degree; index++) {
		sum.coeffs[index] = 0;
		if (index <= P.degree)
			sum.coeffs[index] += P.coeffs[index];
		if (index <= Q.degree)
			sum.coeffs[index] += Q.coeffs[index];
	}
    
	return sum;
}