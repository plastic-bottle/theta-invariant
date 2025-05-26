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
#include <math.h>

 /* The maximum allowable number of crossings a knot can have for the computation of the theta invariant */
#define MAX_CROSSINGS (32)

/* Macros to find min and max of two numbers */
#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

/*Symbolic constants for yes and no*/
enum { YES, NO };

/* Functions to safely use malloc, calloc, and free */
void* safe_malloc(size_t);
void* safe_calloc(size_t, size_t);
void safe_free(void*);

/* Struct for regular polynomial */
/* coeffs[n] stores the coefficient of x^n */
struct polynomial {
	int degree;
	int* coeffs;
};

/* Struct for laurent polynomial */
/* Stores coefficients, highest degree, and lowest degree of polyonmial */
struct laurent_polynomial {
	signed int lowest_degree;
	signed int highest_degree;
	signed int* coeffs;
};

/*Struct for crossing in PD notation; first entry of data is the undercrossing which points at
the crossing (when the knot is given an orientation), and then lists arcs in clockwise order
(this is different from many other PD codes, which go in counterclockwise order) */
struct crossing {
	int data[4];
};

/*Struct for knot in PD notation; contains number of crossings and data of all crossings */
struct knot {
	int number_of_crossings;
	struct crossing* crossings;
};

struct knot make_knot(int, struct crossing*);
int* rotation_numbers(struct knot*);

/*Struct for linked list of integers; stores integer value and pointers to previous/next elements*/
struct linked_list {
	int value;
	struct linked_list* previous;
	struct linked_list* next;
};

struct linked_list* insert_linked_list(int, struct linked_list*, struct linked_list*);
void delete_linked_list(struct linked_list*);
#endif
