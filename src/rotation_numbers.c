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

/*Function which returns an array of rotation numbers for a knot; based off Bar-Natan's code in Mathematica*/
int* rotation_numbers(const struct knot* const K) 
{
	int n = K->number_of_crossings;

	/*Strand positions are 1 and 2 for incoming over/under crossings at a positive crossing, and -1 and -2 for incoming
	over/udner crossings at a negative crossing*/
	int* strand_positions = (int*)safe_malloc((2 * (size_t)n + 1) * sizeof(int));
	/*Strand pair contains the other strand which is also an incoming strand at the current strand's incoming crossing*/
	int* strand_pairs = (int*)safe_malloc((2 * (size_t)n + 1) * sizeof(int));
	for (int crossing = 0; crossing < n; crossing++) {
		int incoming_under = K->crossings[crossing].data[0];
		int incoming_over;
		if ((K->crossings[crossing].data[3] - K->crossings[crossing].data[1] + 2 * n) % (2 * n) == 1) {
			incoming_over = K->crossings[crossing].data[1];
			strand_positions[incoming_over] = 1;
			strand_positions[incoming_under] = 2;
			strand_pairs[incoming_over] = incoming_under;
			strand_pairs[incoming_under] = incoming_over;
		}
		else {
			incoming_over = K->crossings[crossing].data[3];
			strand_positions[incoming_over] = -1;
			strand_positions[incoming_under] = -2;
			strand_pairs[incoming_over] = incoming_under;
			strand_pairs[incoming_under] = incoming_over;
		}
	}

	/*Rotation numbers are calculated as done below*/
	int* rotation_numbers = (int*)safe_calloc((2 * (size_t)n + 1), sizeof(int));
	struct linked_list* front_strand = insert_linked_list(1, NULL, NULL);
	for (int strand = 1; strand <= 2 * n; strand++) {
		int negative_occurs = NO;
		struct linked_list* current_element = front_strand;
		while (current_element != NULL) {
			if (current_element->value == -strand) {
				negative_occurs = YES;
				break;
			}
			current_element = current_element->next;
		}
		current_element = front_strand;
		if (negative_occurs == YES) {
			int last_occurrence = 0;
			while (current_element != NULL) {
				if (current_element->value == -strand && last_occurrence == strand)
					rotation_numbers[strand]--;
				if (abs(current_element->value) == strand)
					last_occurrence = current_element->value;
				current_element = current_element->next;
			}
		}
		else {
			current_element = front_strand;
			while (current_element != NULL) {
				if (current_element->value == strand) {
					delete_linked_list(current_element);
					int pair = strand_pairs[strand];
					struct linked_list* LL;
					if (strand_positions[strand] == 1 || strand_positions[strand] == -2) {
						LL = insert_linked_list(-pair, current_element->previous, current_element->next);
						LL = insert_linked_list(strand + 1, current_element->previous, LL);
						LL = insert_linked_list(pair + 1, current_element->previous, LL);
					}
					else {
						rotation_numbers[pair]++;
						LL = insert_linked_list(pair + 1, current_element->previous, current_element->next);
						LL = insert_linked_list(strand + 1, current_element->previous, LL);
						LL = insert_linked_list(-pair, current_element->previous, LL);
					}
					if (current_element == front_strand)
						front_strand = LL;
				}
				current_element = current_element->next;
			}
		}
	}
	return rotation_numbers;
}
