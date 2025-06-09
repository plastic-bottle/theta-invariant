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

/*Function to print a double polynomial in a certain variable*/
void print_double_polynomial(struct double_polynomial P, char c) {
	for (int degree = P.degree; degree >= 0; degree--) {
		double coeff = P.coeffs[degree];
		if (coeff == 0) {
			if (degree == P.degree)
				printf("0");
			continue;
		}
		if (coeff > 0) {
			if (degree != P.degree)
				printf(" + ");
			if (coeff != 1)
				printf("%f", coeff);
		}
		else {
			printf(" - ");
			if (coeff != -1)
				printf("%f", -coeff);
		}
		if (degree == 0) {
			if (coeff == 1 || coeff == -1)
				printf("1");
		}
		else {
			putchar(c);
			if (degree != 1)
				printf("^%d", degree);
		}
	}
	printf("\n");
}