/* GQSS alignment format function definitions.
 *
 * Copyright (C) 2019 Qijia (Michael) Jin
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef GQSS_ALIGNMENT_FORMAT_H
#define GQSS_ALIGNMENT_FORMAT_H

#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

/*
	generate_int_linear_gap_penalty_pair_alignment(char* substitution_matrix_name, char* query_sequence_identifier, char* sequence_identifier, char* trace_X, char* trace_Y, int64_t score, int64_t gap_penalty)
	
	generate_int_linear_gap_penalty_pair_alignment() returns a formatted pair alignment as a newly allocated C string. The function assumes the alignment's
	linear gap penalty is an integer value.
	
	generate_int_linear_gap_penalty_pair_alignment() will return a NULL pointer if it encounters errors or errorneously formatted function arguments.

	The length of the returned C string was computed using the following numbers:

	Start of Header
	41 + 43 + 37 + 22 + 41 + 41 = (3 x 41) + (43 + 37) + 22 = 123 + 80 + 22 = 123 + 102 = 225
	
	Aligned Sequence Names
	(2 + 23 + 5 + strlen(sequence_id_token + 1) + 1) + (5 + strlen(query_sequence_identifier) + 1)
	= (31 + strlen(sequence_id_token + 1)) + (6 + strlen(query_sequence_identifier))
	= 37 + strlen(sequence_id_token + 1) + strlen(query_sequence_identifier)
	
	Substitution Matrix name
	10 + strlen(substitution_matrix_name) + 1 = 11 + strlen(substitution_matrix_name)
	
	Gap Penalties
	38 + 41 + 2 + 31 = 38 + 43 + 31 = 81 + 31 = 112
	
	Statistics
	65 + 65 + 65 + 65 + 30 = (4 x 65) + 30 = 260 + 30 = 290
	
	End of Header
	2 + 2 + 41 = 4 + 41 = 45
	
	max_seq_id = max(strlen(sequence_id_token + 1), strlen(query_sequence_identifier))

	Section Structure
			2 newline characters			//2
			"%-???s %20llu %s %20llu\n"		//max_seq_id + 1 + 20 + 1 + 50 + 1 + 20 + 1 = max_seq_id + 94
			"%-???s %s\n"					//max_seq_id + 22 + 50 + 1 = max_seq_id + 73
			"%-???s %20llu %s %20llu\n"		//max_seq_id + 1 + 20 + 1 + 50 + 1 + 20 + 1 = max_seq_id + 94
				
		=alignment_quotient x (3 x max() + 2 + 94 + 73 + 94) = alignment_quotient x (3 x max_seq_id + 263)
	
	Remainder Structure
			2 newline characters			//2
			"%-???s %20llu %s %20llu\n"		//max_seq_id + 1 + 20 + 1 + alignment_remaining + 1 + 20 + 1 = max_seq_id + 44 + alignment_remaining
			"%-???s %s\n"					//max_seq_id + 22 + alignment_remaining + 1 = max_seq_id + 23 + alignment_remaining
			"%-???s %20llu %s %20llu\n"		//max_seq_id + 1 + 20 + 1 + alignment_remaining + 1 + 20 + 1 = max_seq_id + 44 + alignment_remaining
				
		=(3 x (max() + alignment_remaining) + 2 + 44 + 23 + 44) = (3 x (max_seq_id + alignment_remaining) + 113)
	
	Footer
	1 + 1 + 41 + 41 = 84
	
	Null Character
	1
*/
char* generate_int_linear_gap_penalty_pair_alignment(char* substitution_matrix_name, char* query_sequence_identifier, char* sequence_identifier, char* trace_X, char* trace_Y, int64_t score, int64_t gap_penalty);

#endif /* GQSS_ALIGNMENT_FORMAT_H */
