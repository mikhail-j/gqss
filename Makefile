# SPDX-License-Identifier: GPL-2.0

CC=gcc

.PHONY: ednafull_linear

ednafull_linear: 
	$(CC) -std=c99 -O2 -o ednafull_linear_smith_waterman linear_gap_smith_waterman.c gqss_file_io.c ednafull_linear_smith_waterman.c

example:
	$(CC) -std=c99 -O2 -o example_linear_gap_smith_waterman linear_gap_smith_waterman.c example_linear_gap_smith_waterman.c
