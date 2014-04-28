/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2014 University of Tennessee / UT-Battelle

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
*******************************************************************************/

#ifndef _SEQ_H
#define _SEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bitmap.h"

#define MAX_SEQ 32000000
#define MAX_LINE 10000
#define WINDOW 120
#define ATG 0
#define GTG 1
#define TTG 2
#define STOP 3
#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000
#define IDEAL_AVG_CONTIG_LEN 1500

int read_seq_training(FILE *, unsigned char *, unsigned char *, double *, int,
                      int *);
int next_seq_multi(FILE *, unsigned char *, unsigned char *, int *, double *,
                   char *, char *);
void rcom_seq(unsigned char *, unsigned char *, unsigned char *, int);

void calc_short_header(char *header, char *short_header, int);

int is_a(unsigned char *, int);
int is_c(unsigned char *, int);
int is_g(unsigned char *, int);
int is_t(unsigned char *, int);
int is_n(unsigned char *, int);
int is_gc(unsigned char *, int);

int is_stop(unsigned char *, int, int);
int is_start(unsigned char *, int, int);
int is_atg(unsigned char *, int);
int is_gtg(unsigned char *, int);
int is_ttg(unsigned char *, int);
int is_nnn(unsigned char *, int);

int codon_has_n(unsigned char *, int);
int gap_to_left(unsigned char *, int);
int gap_to_right(unsigned char *, int);

double prob_stop(int, double);

double gc_content(unsigned char *, int, int);

char amino(unsigned char *, int, int, int);
int amino_num(char);
char amino_letter(int);

int rframe(int, int);
int max_fr(int, int, int);

int *calc_most_gc_frame(unsigned char *, int);

int mer_ndx(int, unsigned char *, int);
void mer_text(char *, int, int);
void calc_mer_bg(int, unsigned char *, unsigned char *, int, double *);

int shine_dalgarno_exact(unsigned char *, int, int, double *);
int shine_dalgarno_mm(unsigned char *, int, int, double *);

int imin(int, int);

#endif
