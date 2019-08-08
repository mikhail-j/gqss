/* Smith-Waterman linear gap penalty example application
 * 
 * This example uses a linear gap penalty of 2 and substitution matrix s(a, b)
 * 
 * where s(a_{i}, b_{j}) = { +3, if a_{i} == b_{j}
 *                         { -3, if a_{i} != b_{j}
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "linear_gap_smith_waterman.h"

#define LINEAR_GAP_PENALTY 2

int64_t get_example_substitution(char a, char b) {
	return (int64_t)(((int)(a == b) * 6) - 3);
}

int main(int argc, char* argv[]) {
	char a[10] = "GGTTGACTA";
	char b[9] = "TGTTACGG";

	//allocate an array of 64-bit integers for our scoring matrix
	int64_t* scores = (int64_t *)malloc(strlen(a) * strlen(b) * sizeof(int64_t));

	//use the Smith-Waterman algorithm
	linear_gap_smith_waterman(a, b, scores, get_example_substitution, LINEAR_GAP_PENALTY);

	//print the result scoring matrix
	printf("Scoring Matrix:\n");
	for (size_t i = 0; i < strlen(a); i++) {
		for (size_t j = 0; j < strlen(b); j++) {
			printf("%2lld ", scores[(i * strlen(b)) + j]);
		}
		printf("\n");
	}

	//obtain the best score of the matrix before tracing our path backwards
	size_t best_i;
	size_t best_j;

	//check if 'best_i' and 'best_j' was set
	//note: please recall that C/C++ is 0-indexed (if using 'best_linear_gap_smith_waterman_score_indices()' from Julia).
	assert(best_linear_gap_smith_waterman_score_indices(strlen(a), strlen(b), scores, &best_i, &best_j));

	//print the matrix indices of the highest score encountered within the matrix
	printf("Best Indices: (%llu, %llu)\n", (uint64_t)best_i, (uint64_t)best_j);

	//allocate C string allocations to store our sequence alignments
	char* trace_a = (char *)malloc((best_i + best_j + 3) * sizeof(char));
	char* trace_b = (char *)malloc((best_i + best_j + 3) * sizeof(char));

	trace_linear_gap_smith_waterman(a, b, scores, trace_a, trace_b, &best_i, &best_j, get_example_substitution, LINEAR_GAP_PENALTY);

	//print in-place updated matrix indices
	printf("Best Indices: (%llu, %llu)\n", (uint64_t)best_i, (uint64_t)best_j);

	//print sequence alignments obtained from our traceback step
	printf("Alignments:\n%s\n%s\n", trace_a, trace_b);

	//free C string allocations
	free(trace_b);
	free(trace_a);

	//free scoring matrix allocation
	free(scores);
	return 0;
}
