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
#include "table.h"

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

double gc_content(unsigned char *, int, int);

int amino_num(char);
char amino_letter(int);

int rframe(int, int);
int max_fr(int, int, int);

int *calc_most_gc_frame(unsigned char *, int);

void mer_text(char *, int, int);
void calc_mer_bg(int, unsigned char *, unsigned char *, int, double *);

int shine_dalgarno_exact(unsigned char *, int, int, double *);
int shine_dalgarno_mm(unsigned char *, int, int, double *);

/*** START INLINE ROUTINES (THOSE TOO SIMPLE TO NOT INLINE) */
/* Simple routines to say whether or not bases are */
/* a, c, t, g, starts, stops, etc. */
inline static int is_a(unsigned char *seq, int n) {
  return nuc(seq, n) == P_A;
}

inline static int is_c(unsigned char *seq, int n) {
  return nuc(seq, n) == P_C;
}

inline static int is_g(unsigned char *seq, int n) {
  return nuc(seq, n) == P_G;
}

inline static int is_t(unsigned char *seq, int n) {
  return nuc(seq, n) == P_T;
}

inline static int is_n(unsigned char *useq, int n) {
  if(test(useq, n) == 0) return 0;
  return 1;
}

inline static int is_stop(unsigned char *seq, int n, struct _training *tinf) {
  unsigned char codon = trinuc(seq, n);
  /* We will worry about "might be a stop or not" at a later date. */
  return (tinf->table[codon] == '*');
}

inline static int is_start(unsigned char *seq, int n, struct _training *tinf) {
  unsigned char codon = trinuc(seq, n);
  return (tinf->table[codon] | 0x20) && (tinf->table[codon] != '*');
}

inline static int is_atg(unsigned char *seq, int n) {
  if(is_a(seq, n) == 0 || is_t(seq, n+1) == 0 || is_g(seq, n+2) == 0) return 0;
  return 1;
}

inline static int is_gtg(unsigned char *seq, int n) {
  if(is_g(seq, n) == 0 || is_t(seq, n+1) == 0 || is_g(seq, n+2) == 0) return 0;
  return 1;
}

inline static int is_ttg(unsigned char *seq, int n) {
  if(is_t(seq, n) == 0 || is_t(seq, n+1) == 0 || is_g(seq, n+2) == 0) return 0;
  return 1;
}

inline static int is_gc(unsigned char *seq, int n) {
  int ndx = n*2;
  if(test(seq, ndx) != test(seq, ndx+1)) return 1;
  return 0;
}

/* Returns a single amino acid for this position */
inline static char amino(unsigned char *seq, int n, struct _training *tinf, int is_init) {
  unsigned char codon = trinuc(seq, n);
  if(is_stop(seq, n, tinf)) return '*';
  if(is_start(seq, n, tinf) && is_init) return 'M';
  return (tinf->table[codon] & (~0xA0));
}

/* Converts a word of size len to a number (opposite bit-order of trinuc!) */
inline static int mer_ndx(int len, unsigned char *seq, int pos) {
  int i, ndx = 0;
  for(i = 0; i < 2*len; i++) ndx |= (test(seq, pos*2+i)<<i);
  return ndx;
}

/* Returns the minimum of two numbers */
inline static int imin(int x, int y) {
  if(x < y) return x;
  return y;
}

#endif
