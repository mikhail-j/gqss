
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
