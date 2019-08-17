/* Smith-Waterman algorithm with a linear gap penalty using the EDNAFULL 
 * substitution matrix. 
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

#include "ednafull_linear_smith_waterman.h"

static const struct option getopt_long_options[] = {
	{"query", required_argument, NULL, 'q'},
	{"gap-penalty", required_argument, NULL, 'P'},
	{"type", required_argument, NULL, 0},
	{"help", no_argument, NULL, 'h'},
	{"version", no_argument, NULL, 'v'},
	{ NULL, 0, NULL, 0}
};

static const char VERSION_STRING[42] = "ednafull_linear_smith_waterman 1.0.0\n";

static const char HELP_STRING[] = (
	"Usage: ednafull_linear_smith_waterman [OPTIONS...] [FASTQ FILE]\n"
	"Run the Smith-Waterman algorithm with linear gap penalty and the EDNAFULL\n"
	"substitution matrix on the given sequences found in the FASTA and FASTQ files.\n\n"
	"Examples:\n"
	"  ednafull_linear_smith_waterman -q gene.fasta reads.fastq\n"
	"  ednafull_linear_smith_waterman -q gene.fasta -P 10 reads.fastq\n"
	"  ednafull_linear_smith_waterman -q gene.fasta --type=pair reads.fastq\n"
	"\n"
	"Options:\n"
	"  -q, --query=FILE            specify query sequence (FASTA format)\n"
	"  -P, --gap-penalty=INT       specify linear gap penalty (default value is 16)\n"
	"  --type=TYPE                 specify output format: 'tsv' (default) or 'pair'\n"
	"  -h, --help                  print this help and exit\n"
	"  --version                   print version information and exit\n"
	);

/* 
	EDNAFULL_NUC_4_4 is actually 90x90 size array.

	The dimensions of array were chosen based on the value of 'Y' (89).
*/
static const int64_t EDNAFULL_NUC_4_4[8100] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, -4, -4, -1, 0, 0, -4, -1, 0, 0, -4, 0, 1, -2, 0, 0, 0, 1, -4, -4, 0, -1, 1, 0, -4,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, -1, -2, 0, 0, -1, -2, 0, 0, -1, 0, -3, -1, 0, 0, 0, -3, -1, -1, 0, -2, -3, 0, -1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, 5, -4, 0, 0, -4, -1, 0, 0, -4, 0, 1, -2, 0, 0, 0, -4, 1, -4, 0, -1, -4, 0, 1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -4, -1, 0, 0, -1, -2, 0, 0, -1, 0, -3, -1, 0, 0, 0, -1, -3, -1, 0, -2, -1, 0, -3,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, -4, -1, 0, 0, 5, -4, 0, 0, 1, 0, -4, -2, 0, 0, 0, 1, 1, -4, 0, -1, -4, 0, -4,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -1, -2, 0, 0, -4, -1, 0, 0, -3, 0, -1, -1, 0, 0, 0, -3, -3, -1, 0, -2, -1, 0, -1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, -4, -1, 0, 0, 1, -3, 0, 0, -1, 0, -4, -1, 0, 0, 0, -2, -2, 1, 0, -3, -2, 0, -2,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -3, 1, -3, 0, 0, -4, -1, 0, 0, -4, 0, -1, -1, 0, 0, 0, -2, -2, -4, 0, -1, -2, 0, -2,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, -2, -1, 0, 0, -1, 0, -1, -1, 0, 0, 0, -1, -1, -2, 0, -1, -1, 0, -1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -3, -4, -1, 0, 0, 1, -3, 0, 0, -2, 0, -2, -1, 0, 0, 0, -1, -2, -4, 0, -1, -2, 0, -4,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, 1, -3, 0, 0, 1, -3, 0, 0, -2, 0, -2, -1, 0, 0, 0, -2, -1, -4, 0, -1, -4, 0, -2,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, -4, -1, 0, 0, -4, -1, 0, 0, 1, 0, -4, -2, 0, 0, 0, -4, -4, 5, 0, -4, 1, 0, 1,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -2, -1, -2, 0, 0, -1, -2, 0, 0, -3, 0, -1, -1, 0, 0, 0, -1, -1, -4, 0, -1, -3, 0, -3,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -3, -4, -1, 0, 0, -4, -1, 0, 0, -2, 0, -2, -1, 0, 0, 0, -2, -4, 1, 0, -3, -1, 0, -2,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -1, 1, -3, 0, 0, -4, -1, 0, 0, -2, 0, -2, -1, 0, 0, 0, -4, -2, 1, 0, -3, -2, 0, -1,
};

/*
	char complement_dna_base(char base)

	complement_dna_base() returns the complement of a given base if possible.
	Otherwise, the function returns the null character.
*/
static char complement_dna_base(char base) {
	switch (base) {
		case 'A':
			return 'T';
		case 'a':
			return 't';
		case 'B':
			return 'V';
		case 'b':
			return 'v';
		case 'C':
			return 'G';
		case 'c':
			return 'g';
		case 'D':
			return 'H';
		case 'd':
			return 'h';
		case 'G':
			return 'C';
		case 'g':
			return 'c';
		case 'H':
			return 'D';
		case 'h':
			return 'd';
		case 'M':
			return 'K';
		case 'm':
			return 'k';
		case 'N':
			return 'N';
		case 'n':
			return 'n';
		case 'S':
			return 'S';
		case 's':
			return 's';
		case 'T':
			return 'A';
		case 't':
			return 'a';
		case 'U':
			return 'A';
		case 'u':
			return 'a';
		case 'V':
			return 'B';
		case 'v':
			return 'b';
		case 'W':
			return 'W';
		case 'w':
			return 'w';
		case 'Y':
			return 'R';
		case 'y':
			return 'r';
		default:
			printf("error: complement_dna_base(): found unexpected base, %c!\n", base);
			return '\0';
	}
}

/*
	char* get_reverse_complement(char* sequence)

	get_reverse_complement() returns the reverse complement of a string in a newly allocated C string.
*/
static char* get_reverse_complement(char* sequence) {
	if (sequence == NULL) {
		return NULL;
	}
	size_t sequence_length = strlen(sequence);

	//allocate a C string with the same size as 'sequence'
	char* reverse_complement_sequence = (char *)malloc((sequence_length + 1) * sizeof(char));
	reverse_complement_sequence[sequence_length] = '\0';

	for (size_t i = 0; i < sequence_length; i++) {
		reverse_complement_sequence[(sequence_length - 1) - i] = complement_dna_base(sequence[i]);
	}

	return reverse_complement_sequence;
}

/*
	int64_t get_nuc_4_4_value(char a, char b)

	int64_t get_nuc_4_4_value(char a, char b) return the value of the 2 bases according to the EDNAFULL substitution matrix.
*/
int64_t get_nuc_4_4_value(char a, char b) {
	size_t index = (size_t)(a) + (90 * (size_t)b);
	return EDNAFULL_NUC_4_4[index];
}

/*
	int64_t get_linear_gap_smith_waterman_score(char* seq_X, char* seq_Y, char** trace_X, char** trace_Y, size_t* start_X, size_t* start_Y, size_t* stop_X, size_t* stop_Y, int64_t gap_penalty)

	get_linear_gap_smith_waterman_score() executes the Smith-Waterman algorithm with linear gap penalty 'gap_penalty' and returns the best score in the matrix.
	The function also sets 'trace_X' and 'trace_Y' to newly allocated C strings that contain the alignment strings. In addition, the indices of the substring are stored into
	'start_X', 'start_Y', 'stop_X', and 'stop_Y'.
*/
int64_t get_linear_gap_smith_waterman_score(char* seq_X, char* seq_Y, char** trace_X, char** trace_Y, size_t* start_X, size_t* start_Y, size_t* stop_X, size_t* stop_Y, int64_t gap_penalty) {
	int64_t score;
	int64_t* Z = (int64_t *)malloc(strlen(seq_X) * strlen(seq_Y) * sizeof(int64_t));

	linear_gap_smith_waterman(seq_X, seq_Y, Z, get_nuc_4_4_value, gap_penalty);

	assert(best_linear_gap_smith_waterman_score_indices(strlen(seq_X), strlen(seq_Y), Z, stop_X, stop_Y));

	//assign initial indices for traceback to 'start_X' and 'start_Y'
	*start_X = *stop_X;
	*start_Y = *stop_Y;

	//allocate alignment strings, this function will not free these allocations
	*trace_X = (char *)malloc(((*stop_X) + (*stop_Y) + 3) * sizeof(char));
	*trace_Y = (char *)malloc(((*stop_X) + (*stop_Y) + 3) * sizeof(char));

	trace_linear_gap_smith_waterman(seq_X, seq_Y, Z, *trace_X, *trace_Y, start_X, start_Y, get_nuc_4_4_value, gap_penalty);

	score = Z[((*stop_X) * strlen(seq_Y)) + *stop_Y];

	//free allocations
	free(Z);

	return score;
}

/*
	double compute_time_elapsed(struct timespec* t1, struct timespec* t2)

	compute_time_elapsed() returns the number of seconds between 2 points in time as a 'double'.
*/
static double compute_time_elapsed(struct timespec* t1, struct timespec* t2) {
	long long int nanoseconds = ((long long int)(t2->tv_sec - t1->tv_sec)*1000000000LL
								+ (long long int)(t2->tv_nsec - t1->tv_nsec));
	return nanoseconds * 0.000000001;
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
	void handle_fastq_tsv(char* fastq_filename, char* fastq_data, char* query_sequence_identifier, char* query_sequence, int64_t gap_penalty)

	handle_fastq_tsv() parses the FASTQ file and writes the results in a tab delimited values file format (TSV).
*/
void handle_fastq_tsv(char* fastq_filename, char* fastq_data, char* query_sequence_identifier, char* query_sequence, int64_t gap_penalty) {
	assert(fastq_filename != NULL);

	size_t total_bytes = strlen(fastq_data);
	size_t current_index = 0;

	uint64_t line_count = 0;
	size_t last_newline = 0;
	size_t current_line_length = 0;

	char* sequence_id = NULL;
	char* sequence = NULL;
	char* phred_scores = NULL;

	char* sequence_alignment;
	char* query_sequence_alignment;
	char* alignment_phred_scores = NULL;
	size_t alignment_phred_scores_length;

	char* reverse_complement_sequence = get_reverse_complement(query_sequence);

	int64_t smith_waterman_score;
	int64_t reverse_complement_smith_waterman_score;

	size_t query_sequence_start;
	size_t query_sequence_stop;
	size_t sequence_start;
	size_t sequence_stop;

	uint64_t identicals;
	uint64_t gaps_X;
	uint64_t gaps_Y;
	uint64_t mismatches;

	//keep track of FASTQ format row as a variable
	size_t sequence_row;

	char* new_filename = (char *)malloc((strlen(fastq_filename) + 8) * sizeof(char));
	if (new_filename == NULL) {
		perror("handle_fastq_tsv(): malloc(): error");

		//immediately exit
		exit(1);
	}

	//determine new .tsv filename from FASTQ file name
	memcpy((new_filename + strlen(fastq_filename)), ".sw.tsv", (8 * sizeof(char)));
	memcpy(new_filename, fastq_filename, (strlen(fastq_filename) * sizeof(char)));

	printf("Writing tab separated values to \"%s\"\n", new_filename);

	FILE* file_fd = fopen(new_filename, "wb");
	if (file_fd == NULL) {
		perror("handle_fastq_tsv(): fopen(): error");

		//immediately exit
		exit(2);
	}

	//free filename string allocation
	free(new_filename);

	//start measuring time between sequences
	struct timespec start_time;
	struct timespec current_time;
	double time_elapsed;

	assert(clock_gettime(CLOCK_MONOTONIC, &start_time) == 0);

	//write the .tsv header (column descriptions) to file
	fprintf(file_fd, "%s", "Reference Sequence Identifier\tSequence Identifier\tSmith-Waterman Score\tLinear Gap Penalty\tSubstitution Matrix\tAlignment Length\tAlignment Identities\tAlignment Gaps\tAlignment Mismatches\tReference Sequence Alignment\tSequence Alignment\tSequence Alignment Base Quality\n");
	if(ferror(file_fd)) {
		perror("handle_fastq_tsv(): fprintf(): error");

		fclose(file_fd);

		//immediately exit
		exit(2);
	}

	while (current_index < total_bytes) {
		if (fastq_data[current_index] == '\n') {
			line_count++;
			current_line_length = current_index - last_newline;
			last_newline = current_index + 1;

			sequence_row = line_count % 4;
			if (sequence_row == 1) {
				//FASTQ sequence identifier
				sequence_id = extract_line(fastq_data, current_index, current_line_length);
			}
			else if (sequence_row == 2) {
				//FASTQ sequence
				sequence = extract_line(fastq_data, current_index, current_line_length);
			}
			else if (sequence_row == 0) {
				//FASTQ quality scores
				phred_scores = extract_line(fastq_data, current_index, current_line_length);

				//run Smith-Waterman algorithm with linear gap
				smith_waterman_score = get_linear_gap_smith_waterman_score(query_sequence, sequence, &sequence_alignment, &query_sequence_alignment, &query_sequence_start, &sequence_start, &query_sequence_stop, &sequence_stop, gap_penalty);

				/*
					Copy the specific section of the FASTQ phred scores corresponding to the alignment.
					
					Note: strlen(alignment_phred_scores) <= strlen(sequence_alignment) due to possible gap insertions in alignment.
				*/
				alignment_phred_scores_length = (sequence_stop - sequence_start) + 1;
				alignment_phred_scores = (char *)malloc((alignment_phred_scores_length + 1) * sizeof(char));
				if (alignment_phred_scores == NULL) {
					perror("handle_fastq_tsv(): malloc(): error");
			
					//immediately exit
					exit(1);
				}
				alignment_phred_scores[alignment_phred_scores_length] = '\0';
			
				memcpy(alignment_phred_scores, (phred_scores + sequence_start), (alignment_phred_scores_length * sizeof(char)));
				

				//count the number of mismatches and gaps found between 'sequence_alignment' and 'query_sequence_alignment'
				count_mismatches(sequence_alignment, query_sequence_alignment, &identicals, &gaps_X, &gaps_Y, &mismatches);

				//format the row output before writing to file
				fprintf(file_fd, "%s\t%s\t%lld\t%lld\t%s\t%llu\t%llu\t%llu\t%llu\t%s\t%s\t%s\n",
								(query_sequence_identifier + 1),
								sequence_id,
								smith_waterman_score,
								gap_penalty,
								"NUC4.4",
								strlen(sequence_alignment),
								identicals,
								(gaps_X + gaps_Y),
								mismatches,
								sequence_alignment,
								query_sequence_alignment,
								alignment_phred_scores);
				if(ferror(file_fd)) {
					perror("handle_fastq_tsv(): fprintf(): error");
			
					fclose(file_fd);
			
					//immediately exit
					exit(2);
				}

				//flush the file stream
				fflush(file_fd);

				//free sequence alignment string allocations
				free(alignment_phred_scores);
				free(sequence_alignment);
				free(query_sequence_alignment);

				//prevent double free() calls by assigning freed memory pointers to NULL
				alignment_phred_scores = NULL;

				//compute the reverse complement sequence alignment
				reverse_complement_smith_waterman_score = get_linear_gap_smith_waterman_score(reverse_complement_sequence, sequence, &sequence_alignment, &query_sequence_alignment, &query_sequence_start, &sequence_start, &query_sequence_stop, &sequence_stop, gap_penalty);

				/*
					Copy the specific section of the FASTQ phred scores corresponding to the alignment.
					
					Note: strlen(alignment_phred_scores) <= strlen(sequence_alignment) due to possible gap insertions in alignment.
				*/
				alignment_phred_scores_length = (sequence_stop - sequence_start) + 1;
				alignment_phred_scores = (char *)malloc((alignment_phred_scores_length + 1) * sizeof(char));
				if (alignment_phred_scores == NULL) {
					perror("handle_fastq_tsv(): malloc(): error");
				
					//immediately exit
					exit(1);
				}
				alignment_phred_scores[alignment_phred_scores_length] = '\0';
				
				memcpy(alignment_phred_scores, (phred_scores + sequence_start), (alignment_phred_scores_length * sizeof(char)));
				
				//count the number of mismatches and gaps found between 'sequence_alignment' and 'query_sequence_alignment'
				count_mismatches(sequence_alignment, query_sequence_alignment, &identicals, &gaps_X, &gaps_Y, &mismatches);

				//format the row output before writing to file
				fprintf(file_fd, "Reverse_Complement_%s\t%s\t%lld\t%lld\t%s\t%llu\t%llu\t%llu\t%llu\t%s\t%s\t%s\n",
								(query_sequence_identifier + 1),
								sequence_id,
								smith_waterman_score,
								gap_penalty,
								"NUC4.4",
								strlen(sequence_alignment),
								identicals,
								(gaps_X + gaps_Y),
								mismatches,
								sequence_alignment,
								query_sequence_alignment,
								alignment_phred_scores);
				if(ferror(file_fd)) {
					perror("handle_fastq_tsv(): fprintf(): error");
			
					fclose(file_fd);
				
					//immediately exit
					exit(2);
				}

				//flush the file stream
				fflush(file_fd);

				//free sequence alignment string allocations
				free(alignment_phred_scores);
				free(sequence_alignment);
				free(query_sequence_alignment);

				//prevent double free() calls by assigning freed memory pointers to NULL
				alignment_phred_scores = NULL;

				//free memory allocations until next sequence
				free(phred_scores);
				free(sequence);
				free(sequence_id);

				//don't call free() twice
				phred_scores = NULL;
				sequence = NULL;
				sequence_id = NULL;

				if (!(line_count & 0x03ff)) {
					//checkpoint after (1024 / 4) = 256 sequences
					assert(clock_gettime(CLOCK_MONOTONIC, &current_time) == 0);
					time_elapsed = compute_time_elapsed(&start_time, &current_time);
		
					printf("[%11.2lf seconds]: %lld sequences parsed\n", time_elapsed, line_count >> 2);
				}
			}
			//else {}//ignore the third line
		}
		current_index++;
	}

	//close file descriptor
	fclose(file_fd);

	//free C string allocations
	free(reverse_complement_sequence);

	//checkpoint after finishing parsing
	assert(clock_gettime(CLOCK_MONOTONIC, &current_time) == 0);
	time_elapsed = compute_time_elapsed(&start_time, &current_time);
	
	printf("[%11.2lf seconds]: %lld sequences parsed\n", time_elapsed, line_count >> 2);

	return;
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
	void handle_fastq_pair(char* fastq_filename, char* fastq_data, char* query_sequence_identifier, char* query_sequence, int64_t gap_penalty)

	handle_fastq_pair() parses the FASTQ file and writes the results in a pair-wise sequence format (pair).
*/
void handle_fastq_pair(char* fastq_filename, char* fastq_data, char* query_sequence_identifier, char* query_sequence, int64_t gap_penalty) {
	assert(fastq_filename != NULL);

	size_t total_bytes = strlen(fastq_data);
	size_t current_index = 0;

	uint64_t line_count = 0;
	size_t last_newline = 0;
	size_t current_line_length = 0;

	char* sequence_id = NULL;
	char* sequence = NULL;
	char* phred_scores = NULL;

	char* sequence_alignment;
	char* query_sequence_alignment;
	char* alignment_phred_scores = NULL;
	size_t alignment_phred_scores_length;

	char* reverse_complement_sequence = get_reverse_complement(query_sequence);

	int64_t smith_waterman_score;
	int64_t reverse_complement_smith_waterman_score;

	size_t query_sequence_start;
	size_t query_sequence_stop;
	size_t sequence_start;
	size_t sequence_stop;

	//keep track of FASTQ format row as a variable
	size_t sequence_row;

	char* new_filename = (char *)malloc((strlen(fastq_filename) + 8) * sizeof(char));
	if (new_filename == NULL) {
		perror("handle_fastq_pair(): malloc(): error");

		//immediately exit
		exit(1);
	}

	//determine new .tsv filename from FASTQ file name
	memcpy((new_filename + strlen(fastq_filename)), ".sw.pair", (9 * sizeof(char)));
	memcpy(new_filename, fastq_filename, (strlen(fastq_filename) * sizeof(char)));

	printf("Writing pair-wise sequence alignments to \"%s\"\n", new_filename);

	FILE* file_fd = fopen(new_filename, "wb");
	if (file_fd == NULL) {
		perror("handle_fastq_pair(): fopen(): error");

		//immediately exit
		exit(2);
	}

	//free filename string allocation
	free(new_filename);

	//start measuring time between sequences
	struct timespec start_time;
	struct timespec current_time;
	double time_elapsed;

	char* alignment_pair = NULL;

	char* query_sequence_id_token = get_first_string_token_space_delimited(query_sequence_identifier);
	assert(query_sequence_id_token != NULL);

	size_t reverse_complement_query_sequence_identifier_length = 19 + strlen(query_sequence_id_token);
	char* reverse_complement_query_sequence_identifier = (char *)malloc((20 + strlen(query_sequence_id_token)) * sizeof(char));
	if (reverse_complement_query_sequence_identifier == NULL) {
		perror("handle_fastq_pair(): malloc(): error");

		//immediately exit
		exit(1);
	}

	memcpy(reverse_complement_query_sequence_identifier, ">Reverse_Complement_", (20 * sizeof(char)));
	memcpy(reverse_complement_query_sequence_identifier + 20, (query_sequence_id_token + 1), ((strlen(query_sequence_id_token) - 1) * sizeof(char)));
	reverse_complement_query_sequence_identifier[reverse_complement_query_sequence_identifier_length] = '\0';

	//free query sequence identifier token string allocation
	free(query_sequence_id_token);

	assert(clock_gettime(CLOCK_MONOTONIC, &start_time) == 0);

	while (current_index < total_bytes) {
		if (fastq_data[current_index] == '\n') {
			line_count++;
			current_line_length = current_index - last_newline;
			last_newline = current_index + 1;

			sequence_row = line_count % 4;
			if (sequence_row == 1) {
				//FASTQ sequence identifier
				sequence_id = extract_line(fastq_data, current_index, current_line_length);
			}
			else if (sequence_row == 2) {
				//FASTQ sequence
				sequence = extract_line(fastq_data, current_index, current_line_length);
			}
			else if (sequence_row == 0) {
				//FASTQ quality scores
				phred_scores = extract_line(fastq_data, current_index, current_line_length);

				//run Smith-Waterman algorithm with linear gap
				smith_waterman_score = get_linear_gap_smith_waterman_score(query_sequence, sequence, &sequence_alignment, &query_sequence_alignment, &query_sequence_start, &sequence_start, &query_sequence_stop, &sequence_stop, gap_penalty);

				//format the sequence alignment output before writing to file
				alignment_pair = generate_int_linear_gap_penalty_pair_alignment("ednafull_linear_smith_waterman", "NUC.4.4", query_sequence_identifier, sequence_id, query_sequence_alignment, sequence_alignment, smith_waterman_score, gap_penalty);

				fprintf(file_fd, "%s", alignment_pair);
				if(ferror(file_fd)) {
					perror("handle_fastq_pair(): fprintf(): error");
			
					fclose(file_fd);
			
					//immediately exit
					exit(2);
				}

				//free pair-wise sequence alignment C string allocation
				free(alignment_pair);

				//flush the file stream
				fflush(file_fd);

				//free sequence alignment string allocations
				free(sequence_alignment);
				free(query_sequence_alignment);

				//prevent double free() calls by assigning freed memory pointers to NULL
				alignment_pair = NULL;
				sequence_alignment = NULL;
				query_sequence_alignment = NULL;

				//compute the reverse complement sequence alignment
				reverse_complement_smith_waterman_score = get_linear_gap_smith_waterman_score(reverse_complement_sequence, sequence, &sequence_alignment, &query_sequence_alignment, &query_sequence_start, &sequence_start, &query_sequence_stop, &sequence_stop, gap_penalty);

				//format the sequence alignment output before writing to file
				alignment_pair = generate_int_linear_gap_penalty_pair_alignment("ednafull_linear_smith_waterman", "NUC.4.4", reverse_complement_query_sequence_identifier, sequence_id, query_sequence_alignment, sequence_alignment, reverse_complement_smith_waterman_score, gap_penalty);

				fprintf(file_fd, "%s", alignment_pair);
				if(ferror(file_fd)) {
					perror("handle_fastq_pair(): fprintf(): error");
			
					fclose(file_fd);
				
					//immediately exit
					exit(2);
				}

				//flush the file stream
				fflush(file_fd);

				//free pair-wise sequence alignment C string allocation
				free(alignment_pair);

				//free sequence alignment string allocations
				free(sequence_alignment);
				free(query_sequence_alignment);

				//prevent double free() calls by assigning freed memory pointers to NULL
				alignment_pair = NULL;
				sequence_alignment = NULL;
				query_sequence_alignment = NULL;

				//free memory allocations until next sequence
				free(phred_scores);
				free(sequence);
				free(sequence_id);

				//don't call free() twice
				phred_scores = NULL;
				sequence = NULL;
				sequence_id = NULL;

				if (!(line_count & 0x03ff)) {
					//checkpoint after (1024 / 4) = 256 sequences
					assert(clock_gettime(CLOCK_MONOTONIC, &current_time) == 0);
					time_elapsed = compute_time_elapsed(&start_time, &current_time);
		
					printf("[%11.2lf seconds]: %lld sequences parsed\n", time_elapsed, line_count >> 2);
				}
			}
			//else {}//ignore the third line
		}
		current_index++;
	}

	//close file descriptor
	fclose(file_fd);

	//free C string allocations
	free(reverse_complement_sequence);
	free(reverse_complement_query_sequence_identifier);

	//checkpoint after finishing parsing
	assert(clock_gettime(CLOCK_MONOTONIC, &current_time) == 0);
	time_elapsed = compute_time_elapsed(&start_time, &current_time);
	
	printf("[%11.2lf seconds]: %lld sequences parsed\n", time_elapsed, line_count >> 2);

	return;
}

/*
	parse_ednafull_linear_smith_waterman_options(int argc, char* argv[], char** query_sequence, char** sequence, int64_t* gap_penalty)

	parse_ednafull_linear_smith_waterman_options() parses the application's given arguments. This function returns 0 when no
	problems were encountered during parsing. Otherwise, parse_ednafull_linear_smith_waterman_options() returns 1 on failure.
*/
static int parse_ednafull_linear_smith_waterman_options(int argc, char* argv[], char** query_sequence, char** sequence, int64_t* gap_penalty, unsigned int* output_flag) {
	int getopt_index = 0;
	int c;

	*query_sequence = NULL;
	*sequence = NULL;

	while ((c = getopt_long(argc, argv, "q:P:hv", getopt_long_options, &getopt_index)) != -1) {
		switch (c) {
			case 0:
				if (strcmp(getopt_long_options[getopt_index].name, "type") == 0) {
					if (strcmp(optarg, "tsv") == 0) {
						*output_flag = OUTPUT_TSV;
					}
					else if (strcmp(optarg, "pair") == 0) {
						*output_flag = OUTPUT_PAIR;
					}
					else {
						printf("ednafull_linear_smith_waterman: option --type: valid types are 'tsv' and 'pair'.\n");
						printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
						return 1;
					}
				}
				break;
			case 'q':
				//check if query file name is an empty string
				if (strlen(optarg) == 0) {
					printf("ednafull_linear_smith_waterman: option -q, --query: FASTA query file name cannot be an empty string.\n");
					printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
					return 1;
				}
				//assign filename
				*query_sequence = optarg;
				break;
			case 'h':
				printf("%s", HELP_STRING);
				exit(0);
				break;
			case 'v':
				printf("%s", VERSION_STRING);
				exit(0);
				break;
			case 'P':
				//assign given gap penalty
				if (sscanf(optarg, "%lld", gap_penalty) == EOF) {
					printf("ednafull_linear_smith_waterman: option -P, --gap-penalty: could not parse the given integer parameter.");
					printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
					return 1;
				}
				break;
			case '?':
				switch (optopt) {
					case 'q':
						printf("ednafull_linear_smith_waterman: option -q, --query: missing FASTA query file name parameter.\n");
						printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
						return 1;
						break;
					case 'P':
						printf("ednafull_linear_smith_waterman: option -P, --gap-penalty: missing gap penalty parameter.\n");
						printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
						return 1;
						break;
					default:
						printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
						return 1;
						break;
				}
				break;
			default:
				printf("ednafull_linear_smith_waterman: unexpected option: %c\n", c);
				printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
				return 2;
				break;
		}
	}

	if (*query_sequence == NULL) {
		printf("ednafull_linear_smith_waterman: expected query sequence file!\n");
		printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
		return 1;
	}
	
	if (argc - optind == 1) {

		if ((strstr(argv[optind], ".fq") == NULL) && (strstr(argv[optind], ".fastq") == NULL)) {
			printf("ednafull_linear_smith_waterman: could not find expected FASTQ file!\n");
			printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
			return 1;
		}
		*sequence = argv[optind];
	}
	else {
		printf("ednafull_linear_smith_waterman: found unexpected number of arguments!\n");
		printf("Try 'ednafull_linear_smith_waterman --help' for more information.\n");
		return 1;
	}

	return 0;
}

int main(int argc, char* argv[]) {
	int64_t gap_penalty = 16;
	char* sequence_filename;
	char* query_sequence_filename;
	unsigned int output_flag;

	int parse_status = parse_ednafull_linear_smith_waterman_options(argc, argv, &query_sequence_filename, &sequence_filename, &gap_penalty, &output_flag);
	
	if (parse_status == 0) {
		char* fasta_sequence_identifier;
		char* query;
		char* fasta_data = read_file(query_sequence_filename);
		size_t fasta_bytes_parsed = extract_fasta_sequence(fasta_data, &fasta_sequence_identifier, &query);
		if (query == NULL) {
			printf("error: failed to read FASTA query sequence!\n");

			free(fasta_data);
			free(fasta_sequence_identifier);
			return 1;
		}

		printf("Query Sequence Identifier: %s\n", (fasta_sequence_identifier + 1));

		char* data = read_file(sequence_filename);
		if (output_flag == OUTPUT_TSV) {
			handle_fastq_tsv(sequence_filename, data, fasta_sequence_identifier, query, gap_penalty);
		}
		else if (output_flag == OUTPUT_PAIR) {
			handle_fastq_pair(sequence_filename, data, fasta_sequence_identifier, query, gap_penalty);
		}
		else {
			printf("error: no output type found!\n");

			//free allocations
			free(data);
			free(query);
			free(fasta_data);
			free(fasta_sequence_identifier);

			return 1;
		}

		//free allocations
		free(data);
		free(query);
		free(fasta_data);
		free(fasta_sequence_identifier);
	}

	return parse_status;
}
