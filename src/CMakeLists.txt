# CMakeList.txt : CMake project for theta_invariant, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("theta_invariant")

# Add source to this project's executable.
add_executable (theta_invariant  "main.c" "theta_implementation.h" "safe_malloc.c" "safe_free.c" "safe_calloc.c" "make_knot.c" "rotation_numbers.c" "insert_linked_list.c" "delete_linked_list.c" "initialize_double_polynomial.c" "make_double_polynomial.c" "add_double_polynomials.c" "multiply_double_polynomials.c" "print_double_polynomial.c" "lagrange_interpolate.c" "adjust_degree.c" "lagrange_product.c" "synthetic_division.c" "int_adjugate.c" "Lagrange Code + Rational Structs (should test for accuracy, but jupyter didn't flag anything).c" "lagrange_polynomial_adjugate.c" "make_float_matrix.c" "make_int_matrix.c" "Quotient of bivariate polynomials (not tested, but should work).c" "row_echelon_form.c")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET theta_invariant PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
