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

 /* Function to scale a bivariate polynomial by a constant factor */
struct bivariate_polynomial scale_bivariate_polynomial(struct bivariate_polynomial P, THETA_INT scale_factor) {
	for (int i = P.lowest_degree_1; i <= P.highest_degree_1; i++) 
		for (int j = P.lowest_degree_2; j <= P.highest_degree_2; j++) 
			MATRIX_ELEMENT(P.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT) *= scale_factor;
}
