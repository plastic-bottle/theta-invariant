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

 /*Function to print a bivariate polynomial in T1 and T2*/
void print_bivariate_polynomial(struct bivariate_polynomial P) {
	enum boolean has_printed = FALSE;
	for (int i = P.highest_degree_1; i >= P.lowest_degree_1; i--) {
		for (int j = P.highest_degree_2; j >= P.lowest_degree_2; j--) {
			THETA_INT coeff = MATRIX_ELEMENT(P.coeffs, i + DEGREE_SHIFT, j + DEGREE_SHIFT);
			if (coeff == 0) {
				if (P.lowest_degree_1 == 0 && P.highest_degree_1 == 0 && P.lowest_degree_2 == 0 && P.highest_degree_2 == 0)
					printf("0");
				continue;
			}
			if (coeff > 0) {
				if (has_printed)
					printf(" + ");
				if (coeff != 1)
					printf("%lld", coeff);
			}
			else {
				printf(" - ");
				if (coeff != -1)
					printf("%lld", -coeff);
			}
			if (i != 0) {
				printf("T1");
				if (i != 1)
					printf("^%d", i);
			}
			if (j != 0) {
				printf("T2");
				if (j != 1)
					printf("^%d", j);
			}
			if (i == 0 && j == 0) {
				if (coeff == 1 || coeff == -1)
					printf("1");
			}
			has_printed = TRUE;
		}
	}
	printf("\n");
}
