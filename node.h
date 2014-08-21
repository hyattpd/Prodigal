/******************************************************************************
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
******************************************************************************/

#ifndef _NODE_H
#define _NODE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "datatypes.h"
#include "sequence.h"

#define STT_NOD 100000
#define MIN_GENE 90
#define MIN_EDGE_GENE 60
#define MAX_SAM_OVLP 60
#define ST_WINDOW 60
#define OPER_DIST 60
#define EDGE_BONUS 0.74
#define EDGE_UPS -1.00
#define META_PEN 7.5
#define MODE_NORM 0
#define MODE_TRN 1
#define MODE_ANON 2

/* Node Types: We use 0 to 9 as starts and 10 to 19 for stops */
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

int add_nodes(unsigned char *, unsigned char *, unsigned char *, int,
              struct _node *, int, int, int);
void zero_nodes(struct _node *, int);
void check_node_allocation(struct _node **, int);
void reset_node_scores(struct _node *, int);
int compare_nodes(const void *, const void *);
int stopcmp_nodes(const void *, const void *);

void record_overlapping_starts(struct _node *, int, double, int);
void frame_score(int *, struct _node *, int);

void score_nodes(unsigned char *, unsigned char *, int, struct _node *, int,
                 struct _training *, int, int);
void calc_coding_score(unsigned char *, unsigned char *, int, struct _node *,
                       int, int, double, double *);
void calc_orf_gc(unsigned char *, struct _node *, int);
void sd_rbs_score(unsigned char *, unsigned char *, int, struct _node *, int,
                  double *);

void find_best_nonsd_motif(struct _training *, unsigned char *, unsigned
                           char *, int, struct _node *, int);
void score_upstream_composition(unsigned char *, int, struct _node *,
                                struct _training *);

double intergenic_mod(struct _node *, struct _node *, double);

int assign_start_value(unsigned char *seq, int);
int assign_stop_value(unsigned char *seq, int);

double dmax(double, double);
double dmin(double, double);

#endif
