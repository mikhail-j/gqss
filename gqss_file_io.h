/* GQSS File I/O related function definitions
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

#ifndef GQSS_FILE_IO_H
#define GQSS_FILE_IO_H

#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include <stdio.h>
#include <assert.h>

#include <sys/stat.h>

//read_file() returns NULL on failure
char* read_file(char* filename);

//extract_line() returns NULL on failure
char* extract_line(char* data, size_t idx, size_t line_length);

//extract first FASTA sequence in 'fasta_data'
char* extract_query_sequence(char* fasta_data, char** fasta_sequence_identifier);

#endif /* GQSS_FILE_IO_H */
