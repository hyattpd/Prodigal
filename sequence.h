/******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2015 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef _SEQ_H
#define _SEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bitmap.h"

#define MAX_SEQ 32000000
#define MAX_LINE 10000
#define WINDOW 120
#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000
#define MIN_CLOSED_CONTIG_LEN 20000
#define SSTRUCT_SIZE 24

#define START 0
#define STOP 1
#define ATG 0
#define GTG 1
#define TTG 2
#define TAA 0
#define TAG 1
#define TGA 2
#define NONST 3   /* Nonstandard start or stop*/
#define EDGE 4    /* Edge start or stop */

int read_seq_training(FILE *, unsigned char *, unsigned char *, double *,
                      int *);
int next_seq_multi(FILE *, unsigned char *, unsigned char *, int *, double *,
                   char *, char *);
void reverse_seq(unsigned char *, unsigned char *, unsigned char *, int);

void calc_short_header(char *header, char *short_header, int);

int is_a(unsigned char *, int);
int is_c(unsigned char *, int);
int is_g(unsigned char *, int);
int is_t(unsigned char *, int);
int is_n(unsigned char *, int);
int is_gc(unsigned char *, int);

int get_stop_type(unsigned char *, int, int);
int get_start_type(unsigned char *, int, int);
int is_atg(unsigned char *, int);
int is_gtg(unsigned char *, int);
int is_ttg(unsigned char *, int);
int is_taa(unsigned char *, int);
int is_tag(unsigned char *, int);
int is_tga(unsigned char *, int);
int is_nnn(unsigned char *, int);

int codon_has_n(unsigned char *, int);
int gap_to_left(unsigned char *, int);
int gap_to_right(unsigned char *, int);

double gc_content(unsigned char *, int, int);

char amino(unsigned char *, int, int, int);
int amino_num(char);
char amino_letter(int);

int assign_start_value(unsigned char *seq, int);
int assign_stop_value(unsigned char *seq, int);

int rframe(int, int);
int max_frame(int, int, int);
int *calc_most_gc_frame(unsigned char *, int);

int mer_index(int, unsigned char *, int);
void mer_text(char *, int, int);
void get_word_counts(int, unsigned char *, unsigned char *, int, double *);

void record_sequence_context(unsigned char *, int, int, int, int *);

int shine_dalgarno_exact(unsigned char *, int, int, double *);
int shine_dalgarno_mismatch(unsigned char *, int, int, double *);

void zero_sequence(unsigned char *, unsigned char *, unsigned char *, int);

int imin(int, int);
int imax(int, int);
double dmin(double, double);
double dmax(double, double);

#endif
