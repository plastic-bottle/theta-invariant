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

/* Allocates memory for a polynomial split matrix */
struct polynomial_split_matrix* make_polynomial_split_matrix(const size_t rows, const size_t cols)
{
    struct polynomial_split_matrix* result = (struct polynomial_split_matrix*) safe_malloc(sizeof(struct polynomial_split_matrix));
    
    result->rows = rows;
    result->cols = cols;

    result->data = (struct polynomial**) safe_malloc(rows * sizeof(struct polynomial*));
    for (int i = 0; i < rows; i++) {
        data[i] = (struct polynomial*) safe_malloc(cols * sizeof(struct polynomial));
    }

    return result;
}
