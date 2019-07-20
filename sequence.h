/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

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
#include <math.h>
#include "bitmap.h"
#include "training.h"
#include "fptr.h"

#define MAX_SEQ 32000000
#define MAX_LINE 10000
#define WINDOW 120
#define MASK_SIZE 50
#define MAX_MASKS 5000
#define ATG 0
#define GTG 1
#define TTG 2
#define STOP 3
#define ACCEPT "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|"

typedef struct _mask {
  int begin;
  int end;
} mask;

int read_seq_training(fptr, unsigned char *, unsigned char *, double *, int,
                      mask *, int *);
int next_seq_multi(fptr, unsigned char *, unsigned char *, int *, double *,
                   int, mask *, int *, char *, char *);
void rcom_seq(unsigned char *, unsigned char *, unsigned char *, int);

void calc_short_header(char *header, char *short_header, int);

int is_a(unsigned char *, int);
int is_c(unsigned char *, int);
int is_g(unsigned char *, int);
int is_t(unsigned char *, int);
int is_n(unsigned char *, int);
int is_gc(unsigned char *, int);

int is_stop(unsigned char *, int, struct _training *);
int is_start(unsigned char *, int, struct _training *);
int is_atg(unsigned char *, int);
int is_gtg(unsigned char *, int);
int is_ttg(unsigned char *, int);

double gc_content(unsigned char *, int, int);

char amino(unsigned char *, int, struct _training *, int);
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
