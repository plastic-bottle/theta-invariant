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

/*Function to initialize double polynomial to the zero polynomial*/
 struct double_polynomial initialize_double_polynomial(void) {
	struct double_polynomial temp;
	temp.degree = 0;
	temp.coeffs = (int*)safe_malloc(MAX_POLYNOMIAL_SIZE * sizeof(double));
	temp.coeffs[0] = 0.0;
	return temp;
}
