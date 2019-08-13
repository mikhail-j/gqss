/* GQSS alignment format functions.
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

#include "gqss_alignment_format.h"

/*
	max_size_t(size_t a, size_t b)

	Return the greatest unsigned integer of the given 'a' and 'b'.
*/
static size_t max_size_t(size_t a, size_t b) {
	if (a < b) {
		return b;
	}
	else {
		return a;
	}
}

/*
	char * get_first_string_token_space_delimited(char* s)

	get_first_string_token_space_delimited() returns a newly allocated C string with the first token
	found by delimiting by the space character (' '). Otherwise, return NULL pointer.
*/
static char * get_first_string_token_space_delimited(char* s) {
	if (s == NULL) {
		return NULL;
	}

	bool found_first_space = false;
	size_t first_space_character;
	for (size_t i = 0; i < strlen(s); i++) {
		if (s[i] == ' ') {
			first_space_character = i;
			found_first_space = true;
			break;
		}
	}
	if (!found_first_space) {
		//return given string 's', no space character was found
		char* token = (char *)malloc((strlen(s) + 1) * sizeof(char));
		if (token == NULL) {
			perror("get_first_string_token_space_delimited(): malloc(): error");

			return NULL;
		}

		token[strlen(s)] = '\0';
		memcpy(token, s, (strlen(s) * sizeof(char)));
		return token;
	}
	else {
		char* token = (char *)malloc((first_space_character + 1) * sizeof(char));
		if (token == NULL) {
			perror("get_first_string_token_space_delimited(): malloc(): error");

			return NULL;
		}
		
		token[first_space_character] = '\0';
		memcpy(token, s, (first_space_character * sizeof(char)));

		return token;
	}
}

/*
	count_mismatches(char* trace_X, char* trace_Y, uint64_t* identical, uint64_t* gaps_X, uint64_t* gaps_Y, uint64_t* mismatches)

	count_mismatches() counts the number of mismatches and gaps found between the 2 given sequences ('trace_X' and 'trace_Y'). The
	resulting number of mismatches found is assigned to 'mismatches'. In addition, the number of gaps corresponding to 'trace_X' and
	'trace_Y' are assigned to 'gaps_X' and 'gaps_Y' respectively.
*/
static void count_mismatches(char* trace_X, char* trace_Y, uint64_t* identical, uint64_t* gaps_X, uint64_t* gaps_Y, uint64_t* mismatches) {
	assert((trace_X != NULL) && (trace_Y != NULL));
	assert(strlen(trace_X) == strlen(trace_Y));

	*identical = 0;
	*gaps_X = 0;
	*gaps_Y = 0;
	*mismatches = 0;

	for (size_t i = 0; i < strlen(trace_X); i++) {
		if (trace_X[i] == trace_Y[i]) {
			if (trace_X[i] == '-') {
				//both bases in 'trace_X' in 'trace_Y' are gaps
				*gaps_X = (*gaps_X) + 1;
				*gaps_Y = (*gaps_Y) + 1;

				*mismatches = (*mismatches) + 1;
			}
			else {
				*identical = (*identical) + 1;
			}
		}
		else {
			if (trace_X[i] == '-') {
				*gaps_X = (*gaps_X) + 1;
			}
			else if (trace_Y[i] == '-') {
				*gaps_Y = (*gaps_Y) + 1;
			}
			
			*mismatches = (*mismatches) + 1;
		}
	}

	return;
}

/*
	generate_int_linear_gap_penalty_pair_alignment(char* substitution_matrix_name, char* query_sequence_identifier, char* sequence_identifier, char* trace_X, char* trace_Y, int64_t score, int64_t gap_penalty)
	
	generate_int_linear_gap_penalty_pair_alignment() returns a formatted pair alignment as a newly allocated C string. The function assumes the alignment's
	linear gap penalty is an integer value.
	
	generate_int_linear_gap_penalty_pair_alignment() will return a NULL pointer if it encounters errors or errorneously formatted function arguments. 
	
	The length of the returned C string was computed using the following numbers:

	Start of Header
	41 + 43 + 37 + 22 + 41 + 41 = (3 x 41) + (43 + 37) + 22 = 123 + 80 + 22 = 123 + 102 = 225
	
	Aligned Sequence Names
	(2 + 23 + 5 + strlen(sequence_id_token + 1) + 1) + (5 + strlen(query_sequence_id_token + 1) + 1)
	= (31 + strlen(sequence_id_token + 1)) + (6 + strlen(query_sequence_id_token + 1))
	= 37 + strlen(sequence_id_token + 1) + strlen(query_sequence_id_token + 1)
	
	Substitution Matrix name
	10 + strlen(substitution_matrix_name) + 1 = 11 + strlen(substitution_matrix_name)
	
	Gap Penalties
	38 + 41 + 2 + 31 = 38 + 43 + 31 = 81 + 31 = 112
	
	Statistics
	65 + 65 + 65 + 65 + 30 = (4 x 65) + 30 = 260 + 30 = 290
	
	End of Header
	2 + 2 + 41 = 4 + 41 = 45
	
	max_seq_id = max(strlen(sequence_id_token + 1), strlen(query_sequence_id_token + 1))

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
char* generate_int_linear_gap_penalty_pair_alignment(char* substitution_matrix_name, char* query_sequence_identifier, char* sequence_identifier, char* trace_X, char* trace_Y, int64_t score, int64_t gap_penalty) {
	assert((trace_X != NULL) && (trace_Y != NULL) && (substitution_matrix_name != NULL));
	assert(strlen(trace_X) == strlen(trace_Y));

	//get the first string token from sequence identifier
	assert(strlen(sequence_identifier) > 1);
	char* sequence_id_token = get_first_string_token_space_delimited(sequence_identifier);
	assert(sequence_id_token != NULL);

	//get the first string token from query sequence identifier
	assert(strlen(query_sequence_identifier) > 1);
	char* query_sequence_id_token = get_first_string_token_space_delimited(query_sequence_identifier);
	assert(query_sequence_id_token != NULL);

	size_t max_sequence_identifier_length = max_size_t(strlen(sequence_id_token + 1), strlen(query_sequence_id_token + 1));

	uint64_t alignment_length = (uint64_t)strlen(trace_X);
	size_t alignment_quotient = alignment_length / 50;
	uint64_t alignment_remaining = alignment_length % 50;

	char* pair_allocation = (char *)malloc(
		(225
		+ 37 + strlen(sequence_id_token + 1) + strlen(query_sequence_id_token + 1)
		+ 11 + strlen(substitution_matrix_name)
		+ 112
		+ 290
		+ 45
		+ (alignment_quotient * ((3 * max_sequence_identifier_length) + 263))
		+ ((3 * (max_sequence_identifier_length + alignment_remaining)) + 113)
		+ 84
		+ 1) * sizeof(char));
	if (pair_allocation == NULL) {
		perror("generate_int_linear_gap_penalty_pair_alignment(): malloc(): error");

		//free C string allocations
		free(query_sequence_id_token);
		free(sequence_id_token);
		return NULL;
	}

	size_t total_bytes_formatted = 0;
	int bytes_written = 0;

	char* time_string = (char *)malloc(25 * sizeof(char));
	if (time_string == NULL) {
		perror("generate_int_linear_gap_penalty_pair_alignment(): malloc(): error");

		//free C string allocations
		free(query_sequence_id_token);
		free(sequence_id_token);
		free(pair_allocation);
		return NULL;
	}

	time_t now = time(NULL);
	struct tm tm_now;

#if defined(__MINGW32__) || defined(__MINGW64__)
	//MinGW does not define localtime_r()
	localtime_s(&tm_now, &now);
#else
	localtime_r(&now, &tm_now);
#endif	/* defined(__MINGW32__) || defined(__MINGW64__) */

	//format time as human-readable C string
	assert(strftime(time_string, 100, "%a %b %d %H:%M:%S %Y", &tm_now) != 0);

	uint64_t identicals;
	uint64_t gaps_X;
	uint64_t gaps_Y;
	uint64_t mismatches;

	//count the number of mismatches and gaps found between 'trace_X' and 'trace_Y'
	count_mismatches(trace_X, trace_Y, &identicals, &gaps_X, &gaps_Y, &mismatches);

	//start of header
	bytes_written = sprintf(pair_allocation,
			"########################################\n"					//41
			"# Program:  ednafull_linear_smith_waterman\n"					//43
			"# Rundate:  %s\n"												//12 + 24 + 1 = 37
			"# Report_file: stdout\n"										//22
			"########################################\n"					//41
			"#=======================================\n", time_string);		//41

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	//sequence identifiers
    bytes_written = sprintf(pair_allocation + total_bytes_formatted,
    		"#\n"															//2
			"# Aligned_sequences: 2\n"										//23
			"# 1: %s\n"														//5 + strlen(sequence_id_token + 1) + 1
			"# 2: %s\n",													//5 + strlen(query_sequence_id_token + 1) + 1
			(sequence_id_token + 1), (query_sequence_id_token + 1));

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	//substitution matrix name
	bytes_written = sprintf(pair_allocation + total_bytes_formatted,
			"# Matrix: %s\n",								//10 + strlen(matrix_name) + 1 = 11 + strlen(substitution_matrix_name)
			substitution_matrix_name);

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	//gap penalties
	bytes_written = sprintf(pair_allocation + total_bytes_formatted,
			"# Gap_penalty: %" PRId64 ".0\n"				//15 + length(repr(typemin(Int64))) + 2 + 1 = 15 + 20 + 2 + 1 = 38
			"# Extend_penalty: %" PRId64 ".0\n"				//18 + length(repr(typemin(Int64))) + 2 + 1 = 18 + 20 + 2 + 1 = 41
			"#\n"											//2
			"# Length: %llu\n",								//10 + length(repr(typemax(UInt64))) + 1 = 10 + 20 + 1 = 31
			gap_penalty, gap_penalty, alignment_length);

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	//statistics
	bytes_written = sprintf(pair_allocation + total_bytes_formatted,
			"# Identity:   %20llu/%llu (%.1f%%)\n"				//14 + 20 + 1 + 20 + 1 + 1 + 6 + 1 + 1 = 14 + 42 + 9 = 65
			"# Similarity: %20llu/%llu (%.1f%%)\n"				//14 + 20 + 1 + 20 + 1 + 1 + 6 + 1 + 1 = 14 + 42 + 9 = 65
			"# Gaps:       %20llu/%llu (%.1f%%)\n"				//14 + 20 + 1 + 20 + 1 + 1 + 6 + 1 + 1 = 14 + 42 + 9 = 65
			"# Mismatchs:  %20llu/%llu (%.1f%%)\n"				//14 + 20 + 1 + 20 + 1 + 1 + 6 + 1 + 1 = 14 + 42 + 9 = 65
			"# Score: %" PRId64 "\n",							//9 + length(typemin(Int64)) + 1 = 9 + 20 + 1 = 30
			identicals, alignment_length, (((double)identicals)/((double)alignment_length) * 100.0),
			identicals, alignment_length, (((double)identicals)/((double)alignment_length) * 100.0),
			(gaps_X + gaps_Y), alignment_length, (((double)(gaps_X + gaps_Y)/(double)(alignment_length)) * 100.0),
			mismatches, alignment_length, (((double)mismatches)/((double)alignment_length) * 100.0),
			score);

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	//end of header
	bytes_written = sprintf(pair_allocation + total_bytes_formatted,
			"#\n"													//2
			"#\n"													//2
			"#=======================================\n");			//41

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	char* sequence_identifier_format_string = (char *)malloc(42 * sizeof(char));
	if (sequence_identifier_format_string == NULL) {
		perror("generate_int_linear_gap_penalty_pair_alignment(): malloc(): error");

		//free C string allocations
		free(query_sequence_id_token);
		free(sequence_id_token);
		free(time_string);
		free(pair_allocation);
		return NULL;
	}

	//format the format string
	assert(sprintf(sequence_identifier_format_string,
		"%%-%llus %%20llu %%s %%20llu\n",							//2 + length(repr(typemax(UInt64))) + 19 = 2 + 20 + 19 = 41
		max_sequence_identifier_length) > 21);

	char* alignment_buffer = (char *)malloc(51 * sizeof(char));
	if (alignment_buffer == NULL) {
		perror("generate_int_linear_gap_penalty_pair_alignment(): malloc(): error");

		//free C string allocations
		free(sequence_identifier_format_string);
		free(query_sequence_id_token);
		free(sequence_id_token);
	    free(time_string);
		free(pair_allocation);
		return NULL;
	}
	alignment_buffer[50] = '\0';

	uint64_t prev_X = 0;
	uint64_t starting_X = 0;
	uint64_t current_X = 0;
	uint64_t prev_Y = 0;
	uint64_t starting_Y = 0;
	uint64_t current_Y = 0;
	for (size_t i = 0; i < alignment_quotient; i++) {
		for (size_t i_i = (i * 50); i_i < ((i + 1) * 50); i_i++) {
			if (trace_X[i_i] != '-') {
				current_X++;
			}
			if (trace_Y[i_i] != '-') {
				current_Y++;
			}
		}

		//do not increment left counter if the section contains zero matches
		if (current_X > prev_X) {
			starting_X = prev_X + 1;
		}
		else {
			starting_X = prev_X;
		}
		if (current_Y > prev_Y) {
			starting_Y = prev_Y + 1;
		}
		else {
			starting_Y = prev_Y;
		}

		//copy section of 'trace_Y' into 'alignment_buffer'
		memcpy(alignment_buffer, (trace_Y + (i * 50)), (50 * sizeof(char)));

		//append 2 newline characters
		pair_allocation[total_bytes_formatted] = '\n';
		pair_allocation[total_bytes_formatted + 1] = '\n';
		pair_allocation[total_bytes_formatted + 2] = '\0';
		total_bytes_formatted = total_bytes_formatted + 2;
		
		//format section of 'trace_Y'
		bytes_written = sprintf(pair_allocation + total_bytes_formatted,
				sequence_identifier_format_string, (sequence_id_token + 1), starting_Y, alignment_buffer, current_Y);
	
		total_bytes_formatted = total_bytes_formatted + bytes_written;
	
		//write white space characters to offset the matches between alignments correctly
		for (size_t j = 0; j < (max_sequence_identifier_length + 22); j++) {
			pair_allocation[total_bytes_formatted + j] = ' ';
		}
		pair_allocation[total_bytes_formatted + (max_sequence_identifier_length + 22)] = '\0';
	
		total_bytes_formatted = total_bytes_formatted + (max_sequence_identifier_length + 22);

		//indicate where matches occur and end with newline character
		for (size_t j = 0; j < 50; j++) {
			if (trace_X[(i * 50) + j] == trace_Y[(i * 50) + j]) {
				if (trace_X[(i * 50) + j] == '-') {
					pair_allocation[total_bytes_formatted + j] = ' ';
				}
				else {
					pair_allocation[total_bytes_formatted + j] = '|';
				}
			}
			else {
				pair_allocation[total_bytes_formatted + j] = ' ';
			}
		}
		pair_allocation[total_bytes_formatted + 50] = '\n';
		pair_allocation[total_bytes_formatted + 51] = '\0';
	
		total_bytes_formatted = total_bytes_formatted + 51;

		//copy section of 'trace_X' into 'alignment_buffer'
		memcpy(alignment_buffer, (trace_X + (i * 50)), (50 * sizeof(char)));

		//format section of 'trace_X'
		bytes_written = sprintf(pair_allocation + total_bytes_formatted,
			sequence_identifier_format_string, (query_sequence_id_token + 1), starting_X, alignment_buffer, current_X);
	
		total_bytes_formatted = total_bytes_formatted + bytes_written;

		prev_X = current_X;
		prev_Y = current_Y;
	}

	if (alignment_remaining != 0) {
		for (size_t i = (alignment_quotient * 50); i < (alignment_quotient * 50) + alignment_remaining; i++) {
			if (trace_X[i] != '-') {
				current_X++;
			}
			if (trace_Y[i] != '-') {
				current_Y++;
			}
		}

		//do not increment left counter if the section contains zero matches
		if (current_X > prev_X) {
			starting_X = prev_X + 1;
		}
		else {
			starting_X = prev_X;
		}
		if (current_Y > prev_Y) {
			starting_Y = prev_Y + 1;
		}
		else {
			starting_Y = prev_Y;
		}

		alignment_buffer[alignment_remaining] = '\0';

		//copy section of 'trace_Y' into 'alignment_buffer'
		memcpy(alignment_buffer, (trace_Y + (alignment_quotient * 50)), (alignment_remaining * sizeof(char)));
		

		//append 2 newline characters
		pair_allocation[total_bytes_formatted] = '\n';
		pair_allocation[total_bytes_formatted + 1] = '\n';
		pair_allocation[total_bytes_formatted + 2] = '\0';
		total_bytes_formatted = total_bytes_formatted + 2;
		
		//format section of 'trace_Y'
		bytes_written = sprintf(pair_allocation + total_bytes_formatted,
				sequence_identifier_format_string, (sequence_id_token + 1), starting_Y, alignment_buffer, current_Y);

		total_bytes_formatted = total_bytes_formatted + bytes_written;

		//write white space characters to offset the matches between alignments correctly
		for (size_t j = 0; j < (max_sequence_identifier_length + 22); j++) {
			pair_allocation[total_bytes_formatted + j] = ' ';
		}
		pair_allocation[total_bytes_formatted + (max_sequence_identifier_length + 22)] = '\0';
	
		total_bytes_formatted = total_bytes_formatted + (max_sequence_identifier_length + 22);
	
		//indicate where matches occur and end with newline character
			for (size_t j = 0; j < alignment_remaining; j++) {
				if (trace_X[(alignment_quotient * 50) + j] == trace_Y[(alignment_quotient * 50) + j]) {
					if (trace_X[(alignment_quotient * 50) + j] == '-') {
						pair_allocation[total_bytes_formatted + j] = ' ';
					}
					else {
						pair_allocation[total_bytes_formatted + j] = '|';
					}
				}
				else {
					pair_allocation[total_bytes_formatted + j] = ' ';
				}
			}
		pair_allocation[total_bytes_formatted + alignment_remaining] = '\n';
		pair_allocation[total_bytes_formatted + alignment_remaining + 1] = '\0';
	
		total_bytes_formatted = total_bytes_formatted + (alignment_remaining + 1);
		
		//copy section of 'trace_X' into 'alignment_buffer'
		memcpy(alignment_buffer, (trace_X + (alignment_quotient * 50)), (alignment_remaining * sizeof(char)));
		
		//format section of 'trace_X'
		bytes_written = sprintf(pair_allocation + total_bytes_formatted,
				sequence_identifier_format_string, (query_sequence_id_token + 1), starting_X, alignment_buffer, current_X);

		total_bytes_formatted = total_bytes_formatted + bytes_written;
	}

	//footer
	bytes_written = sprintf(pair_allocation + total_bytes_formatted,
			"\n"												//1
			"\n"												//1
			"#---------------------------------------\n"		//41
			"#---------------------------------------\n");		//41

	total_bytes_formatted = total_bytes_formatted + bytes_written;

	//check if the total number of characters matches the C string length
	assert(total_bytes_formatted == strlen(pair_allocation));
	
	//free C string allocations
	free(alignment_buffer);
	free(sequence_identifier_format_string);
	free(query_sequence_id_token);
	free(sequence_id_token);
    free(time_string);
    
	return pair_allocation;
}
