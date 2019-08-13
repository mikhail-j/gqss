/* GQSS File I/O related functions.
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

#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include <stdio.h>
#include <assert.h>

#include <sys/stat.h>

//read_file() returns NULL on failure
char* read_file(char* filename) {
	struct stat file_stat;
	size_t bytes_read = 0;
	size_t fread_bytes;

	//get file size
	assert((stat(filename, &file_stat) == 0) && ((file_stat.st_mode & S_IFMT) == S_IFREG));

	char* file_data = (char *)malloc((file_stat.st_size + 1) * sizeof(char));
	if (file_data == NULL) {
		perror("error: malloc(): memory could not be allocated for file");
		return file_data;
	}
	file_data[file_stat.st_size] = '\0';

	FILE* file_fd = fopen(filename, "rb");
	if (file_fd == NULL) {
		perror("error: fopen()");

		free(file_data);
		file_data = NULL;
		return file_data;
	}

	while (!feof(file_fd)) {
		fread_bytes = fread(file_data + bytes_read, sizeof(char), 1024, file_fd);
		bytes_read = bytes_read + fread_bytes;
	}

	//make sure we read the correct amount of data
	assert(bytes_read == file_stat.st_size);

	fclose(file_fd);
	return file_data;
}

//extract_line() returns NULL on failure
char* extract_line(char* data, size_t idx, size_t line_length) {
	char* line = (char *)malloc((line_length + 1) * sizeof(char));
	if (line == NULL) {
		perror("error: extract_line(): malloc()");
		return line;
	}

	//set null terminator
	line[line_length] = '\0';

	//copy line to new C string allocation
	memcpy(line, (data + (idx - line_length)), line_length * sizeof(char));

	//check for carriage return
	if (line[line_length - 1] == '\r') {
		line[line_length - 1] = '\0';
	}

	return line;
}

//extract first FASTA sequence in 'fasta_data'
//set 'fasta_sequence_identifier' to the sequence identifier corresponding to the 'sequence' returned
char* extract_query_sequence(char* fasta_data, char** fasta_sequence_identifier) {
	char* sequence = NULL;
	*fasta_sequence_identifier = NULL;

	size_t total_bytes = strlen(fasta_data);
	size_t current_index = 0;

	size_t line_count = 0;
	size_t last_newline = 0;
	size_t current_line_length = 0;

	while (current_index < total_bytes) {
		if (fasta_data[current_index] == '\n') {
			line_count++;
			current_line_length = current_index - last_newline;
			last_newline = current_index + 1;

			if (line_count == 1) {
				*fasta_sequence_identifier = extract_line(fasta_data, current_index, current_line_length);
			}
			else if (line_count == 2) {
				sequence = extract_line(fasta_data, current_index, current_line_length);
				break;
			}
		}
		current_index++;
	}

	return sequence;
}
