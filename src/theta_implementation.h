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
/* The maximum number of coefficients in any polynomial used to compute the theta invariant*/
#define MAX_POLYNOMIAL_SIZE (2 * MAX_CROSSINGS + 1)

/* Macros to find min and max of two numbers */
#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

/* Symbolic constants for yes and no */
enum { YES, NO };

/* booleans */
typedef int bool;
enum { false, true };

/* Macro for floating point data type we are using */
#define THETA_FLOAT double
#define THETA_INT long long


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

/* Struct for polynomial with floating-point coefficients */
/* coeffs[n] stores the coefficient of x^n */
struct double_polynomial {
	size_t degree;
	double* coeffs;
};

extern struct polynomial initialize_polynomial(void);

extern struct double_polynomial initialize_double_polynomial(void);
extern struct double_polynomial make_double_polynomial(size_t degree, double* coeffs);
extern struct double_polynomial add_double_polynomials(struct double_polynomial P, struct double_polynomial Q);
extern struct double_polynomial multiply_double_polynomials(struct double_polynomial P, struct double_polynomial Q);
extern void print_double_polynomial(struct double_polynomial P, char c);
extern struct double_polynomial lagrange_product(int max_degree, double* inputs);
extern struct double_polynomial synthetic_division(struct double_polynomial P, double a);
extern struct double_polynomial lagrange_interpolate(int max_degree, double* inputs, double* outputs);
extern void adjust_degree(struct double_polynomial* P);

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

/* Matrix of ints */
struct int_matrix {
	size_t rows;
	size_t cols;
	THETA_INT* data;
};

/* Allocates memory for a polynomial matrix */
extern struct polynomial_matrix* make_polynomial_matrix(const size_t rows, const size_t cols);

/* Allocates memory for a float matrix */
extern struct float_matrix* make_float_matrix(const size_t rows, const size_t cols);

/* Allocates memory for a int matrix */
extern struct int_matrix* make_int_matrix(const size_t rows, const size_t cols);


/* Puts results in order_basis and s_col_degs */
/* If the input and result arguments are identical, then old memory will be freed and new memory will be allocated */
extern void right_order_basis(struct polynomial_matrix* const F, int order, int* const s_array, struct polynomial_matrix* const order_basis, int* const s_col_degs)

/* Helper function for the minimal_nullspace_basis algorithm */
/* We require A and B to satisfy the constraints of Theorem 3.7 in Zhou et al. (2012) */
extern struct polynomial_matrix* mnb_fast_multiplication(struct polynomial_matrix* const A, struct polynomial_matrix* const B);

#define ORDER_CONST 3

/* Requires F to be s-reduced */
/* Puts results in nullspace_basis and s_col_degs */
/* If any input and result arguments are identical, then old memory will be freed and new memory will be allocated */
extern void minimal_nullspace_basis(struct polynomial_matrix* const F, int* const s_array, struct polynomial_matrix* const nullspace_basis, int* const s_col_degs)


/* Puts float matrix A into REF and stores in result without scaling any rows to preserve determinant */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions */
extern void row_echelon_form(struct float_matrix* const A, struct float_matrix* const result);

/* Returns the determinant of A, and populates matrix result with the entries of the adjugate of A */
/* A and result are allowed to point to the same memory */
/* Assumes that A and result have the same dimensions and are both square */
extern THETA_INT int_adjugate(struct int_matrix* const A, struct int_matrix* const result);

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
