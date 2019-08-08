/* Functions that implement the Smith-Waterman algorithm with a linear gap penalty.
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

#include "linear_gap_smith_waterman.h"

/*
	max(int64_t a, int64_t b)

	Return the greatest integer of the given 'a' and 'b'.
*/
static int64_t max(int64_t a, int64_t b) {
	if (a < b) {
		return b;
	}
	else {
		return a;
	}
}

/*
	best_linear_gap_smith_waterman_score(int64_t left, int64_t up_left, int64_t up, char a, char b, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty)

	best_linear_gap_smith_waterman_score() returns the best possible score for an element of the scoring matrix based on the
	neighbors of that element and the characters 'a' and 'b' of 2 sequences to be aligned.
*/
int64_t best_linear_gap_smith_waterman_score(int64_t left, int64_t up_left, int64_t up, char a, char b, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty) {
	return max(max(max(left - gap_penalty, up - gap_penalty), (up_left + get_substitution_matrix_value(a, b))), 0);
}

/*
	linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* scores, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty)

	linear_gap_smith_waterman() is an implementation of the Smith-Waterman algorithm with a linear gap penalty.
*/
void linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* scores, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty) {
	size_t len_X = strlen(seq_X);
	size_t len_Y = strlen(seq_Y);

	//first row done without loop
	scores[0] = best_linear_gap_smith_waterman_score(0, 0, 0, seq_X[0], seq_Y[0], get_substitution_matrix_value, gap_penalty);
	for (size_t j = 1; j < len_Y; j++) {
		scores[j] = best_linear_gap_smith_waterman_score(scores[j - 1], 0, 0, seq_X[0], seq_Y[j], get_substitution_matrix_value, gap_penalty);
	}

	for (size_t i = 1; i < len_X; i++) {
		scores[(i * len_Y)] = best_linear_gap_smith_waterman_score(0,
										0,
										scores[((i - 1) * len_Y)],
										seq_X[i],
										seq_Y[0], get_substitution_matrix_value, gap_penalty);
		for (size_t j = 1; j < len_Y; j++) {
			scores[(i * len_Y) + j] = best_linear_gap_smith_waterman_score(scores[(i * len_Y) + j - 1],
												scores[((i - 1) * len_Y) + j - 1],
												scores[((i - 1) * len_Y) + j],
												seq_X[i],
												seq_Y[j], get_substitution_matrix_value, gap_penalty);
		}
	}
	return;
}

/*
	best_linear_gap_smith_waterman_score_indices(size_t len_X, size_t len_Y, int64_t* Z, size_t* x, size_t* y)

	best_linear_gap_smith_waterman_score_indices() returns true if 'x' and 'y' were assigned the indices of the best
	score in the given matrix Z. Otherwise, return false if the matrix contains no elements.
*/
bool best_linear_gap_smith_waterman_score_indices(size_t len_X, size_t len_Y, int64_t* Z, size_t* x, size_t* y) {
	/*
		Initialize best score to -1 (which is an impossible value of an element for the scoring
		matrix of the Smith-Waterman algorithm).
	*/
	int64_t best_score = -1;

	//Check if matrix is empty
	if ((len_X == 0) || (len_Y == 0)) {
		return false;
	}

	size_t array_index = 0;
	for (size_t i = 0; i < len_X; i++) {
		for (size_t j = 0; j < len_Y; j++) {
			if (Z[array_index] > best_score) {
				best_score = Z[array_index];
				*x = i;
				*y = j;
			}		
			++array_index;
		}
	}

	return true;
}

/*
	trace_linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* Z, char* trace_X, char* trace_Y, size_t x, size_t y, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty)

	trace_linear_gap_smith_waterman() expects a matrix scored by the Smith-Waterman algorithm. The strings 'trace_X' and 'trace_Y'
	should be given alignment 'char *' allocations of size (length(X) + length(Y) + 1) for worst case (triangle inequality)
*/
void trace_linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* Z, char* trace_X, char* trace_Y, size_t* x, size_t* y, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty) {
	size_t len_X = strlen(seq_X);
	size_t len_Y = strlen(seq_Y);
	assert(((len_X > 0) && (len_Y > 0)));


	int64_t score = Z[((*x) * len_Y) + (*y)];

	size_t alignment_index = 0;

	//we should break when we see the next match is 0
	while (score != 0) {
		if ((*x == 0) || (*y == 0)) {
			trace_X[alignment_index] = seq_X[*x];
			trace_Y[alignment_index] = seq_Y[*y];
			break;
		}

		//check left, top/left, top cells
		if (Z[((*x * len_Y) + *y - 1)] - gap_penalty == Z[((*x) * len_Y) + (*y)]) {
			trace_X[alignment_index] = '-';
			trace_Y[alignment_index] = seq_Y[*y];

			score = Z[(((*x) * len_Y) + (*y) - 1)];

			*y = *y - 1;
			alignment_index++;
		}
		else if (Z[((((*x) - 1) * len_Y) + ((*y) - 1))] + get_substitution_matrix_value(seq_X[*x], seq_Y[*y]) == Z[((*x) * len_Y) + (*y)]) {
			trace_X[alignment_index] = seq_X[*x];
			trace_Y[alignment_index] = seq_Y[*y];

			//check if next diagonal cell is zero
			if (Z[(((*x) - 1) * len_Y) + ((*y) - 1)] == 0) {
				break;
			}

			score = Z[((((*x) - 1) * len_Y) + ((*y) - 1))];

			*x = *x - 1;
			*y = *y - 1;
			alignment_index++;
		}
		else if (Z[((((*x) - 1) * len_Y) + *y)] - gap_penalty == Z[((*x) * len_Y) + *y]) {
			trace_X[alignment_index] = seq_X[*x];
			trace_Y[alignment_index] = '-';

			score = Z[((((*x) - 1) * len_Y) + (*y))];

			*x = *x - 1;
			alignment_index++;
		}
		else {
			//we shouldn't reach here!
			assert(false);
		}
	}

	size_t alignment_length = alignment_index + 1;

	trace_X[alignment_length] = '\0';
	trace_Y[alignment_length] = '\0';

	char swap_buffer;
	for (size_t i = 0; i < (alignment_length >> 1); i++) {
		swap_buffer = trace_X[i];
		trace_X[i] = trace_X[alignment_index - i];
		trace_X[alignment_index - i] = swap_buffer;

		swap_buffer = trace_Y[i];
		trace_Y[i] = trace_Y[alignment_index - i];
		trace_Y[alignment_index - i] = swap_buffer;
	}
	return;
}