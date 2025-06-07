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

/* Symbolic constants for yes and no */
enum { YES, NO };

/* Macro for floating point data type we are using */
#define THETA_FLOAT __float128



/* Functions to safely use malloc, calloc, and free */
extern void* safe_malloc(const size_t size);
extern void* safe_calloc(const size_t n, const size_t size);
extern void safe_free(void* allocated_pointer);



/* Struct for regular polynomial */
/* coeffs[n] stores the coefficient of x^n */
struct polynomial {
	size_t degree;
	int* coeffs;
};

/* Struct for laurent polynomial */
/* Stores coefficients, highest degree, and lowest degree of polyonmial */
struct laurent_polynomial {
	signed int lowest_degree;
	signed int highest_degree;
	signed int* coeffs;
};



/* Macro for accessing elements of a pointer to a matrix */
#define MATRIX_ELEMENT(A, row, col) A->data[(size_t) (row) * A->cols + (size_t) (col)]

/* Matrix of polynomial structs */
struct polynomial_matrix {
	size_t rows;
	size_t cols;
	struct polynomial* data;
};

/* Matrix of floats */
struct float_matrix {
	size_t rows;
	size_t cols;
	THETA_FLOAT* data;
};

/* Allocates memory for a float matrix */
extern struct float_matrix* make_float_matrix(const size_t rows, const size_t cols);

/* Puts float matrix A into REF and stores in result without scaling any rows to preserve determinant */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions */
extern void row_echelon_form(struct float_matrix* const A, struct float_matrix* const result);

/* Returns the determinant of A, and populates matrix result with the entries of the adjugate of A */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions and are both square */
extern THETA_FLOAT float_adjugate(struct float_matrix* const A, struct float_matrix* const result);

/* Returns the adjugate of polynomial matrix A via Lagrange interpolation at consecutive integers starting from start_t */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions and are both square */
extern void lagrange_polynomial_adjugate(struct polynomial_matrix* const A, struct polynomial_matrix* const result, const int start_t);



/* Struct for crossing in PD notation; first entry of data is the undercrossing which points at
   the crossing (when the knot is given an orientation), and then lists arcs in clockwise order
   (this is different from many other PD codes, which go in counterclockwise order) */
struct crossing {
	int data[4];
};

/* Struct for knot in PD notation; contains number of crossings and data of all crossings */
struct knot {
	size_t number_of_crossings;
	struct crossing* crossings;
};

extern struct knot make_knot(const int number_of_crossings, struct crossing* const crossings);
extern int* rotation_numbers(const struct knot* const K);



/* Struct for linked list of integers; stores integer value and pointers to previous/next elements */
struct linked_list {
	int value;
	struct linked_list* previous;
	struct linked_list* next;
};

extern struct linked_list* insert_linked_list(int value, struct linked_list* previous, struct linked_list* next);
extern void delete_linked_list(struct linked_list* LL);

#endif
