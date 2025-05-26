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

#ifndef THETA_IMPL_H
#define THETA_IMPL_H

#include <stdio.h>
#include <stdlib.h>

/* The maximum allowable number of crossings a knot can have for the computation of the theta invariant */
#define MAX_CROSSINGS (32)

/* Macros to find min and max of two numbers */
#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

/* Functions to safely use malloc, calloc, and free */
void* safe_malloc(size_t);
void* safe_calloc(size_t, size_t);
void safe_free(void*);

/* Struct for regular polynomial */
/* coeffs[n] stores the coefficient of x^n */
struct polynomial{
    int degree;
    int *coeffs;
};

/* Struct for laurent polynomial */
/* Stores coefficients, highest degree, and lowest degree of polyonmial */
struct laurent_polynomial {
	signed int lowest_degree;
	signed int highest_degree;
	signed int* coeffs;
};


#endif