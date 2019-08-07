/* Function definitions for the Smith-Waterman algorithm with a linear gap penalty.
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

#ifndef GQSS_LINEAR_GAP_SMITH_WATERMAN_H
#define GQSS_LINEAR_GAP_SMITH_WATERMAN_H

#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

/*
	best_linear_gap_smith_waterman_score(int64_t left, int64_t up_left, int64_t up, char a, char b, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty)

	best_linear_gap_smith_waterman_score() returns the best possible score for an element of the scoring matrix based on the
	neighbors of that element and the characters 'a' and 'b' of 2 sequences to be aligned.
*/
int64_t best_linear_gap_smith_waterman_score(int64_t left, int64_t up_left, int64_t up, char a, char b, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty);

/*
	linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* scores, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty)

	linear_gap_smith_waterman() is an implementation of the Smith-Waterman algorithm with a linear gap penalty.
*/
void linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* scores, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty);

/*
	best_linear_gap_smith_waterman_score_indices(size_t len_X, size_t len_Y, int64_t* Z, size_t* x, size_t* y)

	best_linear_gap_smith_waterman_score_indices() returns true if 'x' and 'y' were assigned the indices of the best
	score in the given matrix Z. Otherwise, return false if the matrix contains no elements.
*/
bool best_linear_gap_smith_waterman_score_indices(size_t len_X, size_t len_Y, int64_t* Z, size_t* x, size_t* y);

/*
	trace_linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* Z, char* trace_X, char* trace_Y, size_t x, size_t y, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty)

	trace_linear_gap_smith_waterman() expects a matrix scored by the Smith-Waterman algorithm. The strings 'trace_X' and 'trace_Y'
	should be given alignment 'char *' allocations of size (length(X) + length(Y) + 1) for worst case (triangle inequality)
*/
void trace_linear_gap_smith_waterman(char* seq_X, char* seq_Y, int64_t* Z, char* trace_X, char* trace_Y, size_t x, size_t y, int64_t (*get_substitution_matrix_value)(char a, char b), int64_t gap_penalty);

#endif /* GQSS_LINEAR_GAP_SMITH_WATERMAN_H */
