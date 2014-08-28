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

#ifndef _TRAIN_H
#define _TRAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datatypes.h"
#include "dprog.h"
#include "node.h"
#include "sequence.h"
#include "summary.h"

#define COD_SKIP 24
#define SD_ITER 10
#define NONSD_ITER 20
#define MIN_TRAIN_GENE_SCORE 35.0
#define MIN_AVG_TRAIN_GENE_LEN 600
#define MIN_AVG_TRAIN_CTG_LEN 3000
#define SMALL_TRAIN_SET_NODES 2000

extern void log_text(int, char *);
int write_training_file(char *, struct _training *);
int read_training_file(char *, struct _training *);

void build_training_set_full(struct _node *, struct _training *,
                             struct _summary *, unsigned char *,
                             unsigned char *, unsigned char *, int,
                             int *, int, int, int);
void build_training_set(struct _node *, struct _training *, struct _summary *,
                        unsigned char *, unsigned char *, unsigned char *, int,
                        int *, int);
void record_gc_frame_bias(struct _training *, unsigned char *, int, struct
                          _node *, int);
void calc_dicodon_gene(struct _training *, unsigned char *, unsigned char *,
                       int, struct _node *, int);
int training_set_quality(struct _summary *);
void low_gene_len_warning(int, struct _summary *);

void train_starts_sd(unsigned char *, unsigned char *, int, struct _node *,
                        int, struct _training *);
void train_starts_nonsd(unsigned char *, unsigned char *, int, struct _node *,
                        int, struct _training *);

void label_good_nonsd_motifs(double [4][4][4096], int [4][4][4096], double);
void update_nonsd_motif_counts(double [4][4][4096], double *, unsigned char *,
                               unsigned char *, int, struct _node *, int);

void determine_sd_usage(struct _training *);
void zero_start_weights(struct _training *, int);
void calc_type_background(struct _node *, int, double *);
void calc_dimer_background(struct _node *, int, double *);
void calc_pair_comp_background(struct _node *, int, double *);
void calc_sd_rbs_background(struct _node *, int, double *, double *);

void normalize_array(double *, int);
void create_log_score(double *, double *, double *, int);

#endif
