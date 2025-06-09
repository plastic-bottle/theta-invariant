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

/*Function which returns the unique polynomial of minimal degree passing through a set of points*/
struct double_polynomial lagrange_interpolate(int max_degree, double* inputs, double* outputs) {
    struct double_polynomial result = initialize_double_polynomial();
    for (int index = 0; index <= max_degree; index++) {
        /*Method for computing basis polynomial is inefficient and could be sped up using Evan's
        code for polynomial quotients*/
        struct double_polynomial basis_polynomial = make_double_polynomial(0, outputs + index);
        for (int second_index = 0; second_index <= max_degree; second_index++) {
            if (second_index == index)
                continue;
            double coeffs[] = { -inputs[second_index] / (inputs[index] - inputs[second_index]), 1.0 / (inputs[index] - inputs[second_index]) };
            basis_polynomial = multiply_double_polynomials(basis_polynomial, make_double_polynomial(1, coeffs));
        }
        result = add_double_polynomials(result, basis_polynomial);
    }
    return result;
}