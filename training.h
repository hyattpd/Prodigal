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

#ifndef _TRAIN_H
#define _TRAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


struct _training {
  double gc;                    /* GC Content */
  int trans_table;              /* 11 = Standard Microbial, NCBI Trans Table to
                                   use */
  double st_wt;                 /* Start weight */
  double bias[3];               /* GC frame bias for each of the 3 positions */
  double type_wt[3];            /* Weights for ATG vs GTG vs TTG */
  int uses_sd;                  /* 0 if doesn't use SD motif, 1 if it does */
  double rbs_wt[28];            /* Set of weights for RBS scores */
  double ups_comp[32][4];       /* Base composition weights for non-RBS-distance
                                   motifs.  0-1 are the -1/-2 position, 2-31 are
                                   the -15 to -44 positions.  Second array is
                                   the base A,C,T,G,etc. */
  double mot_wt[4][4][4096];    /* Weights for upstream motifs.  First index is
                                   the motif length (3-6), the second is the
                                   spacer distance (0 = 5-10bp, 1 = 3-4bp, 2 =
                                   11-12bp, 3 = 13-15bp), and the last is the
                                   numerical value of the motif (ranging from 0
                                   to 4095 for 6-mers, less for shorter
                                   motifs) */
  double no_mot;                /* Weight for the case of no motif */
  double gene_dc[4096];         /* Coding statistics for the genome */
};

int write_training_file(char *, struct _training *);
int read_training_file(char *, struct _training *);

void initialize_metagenome_0(struct _training *);
void initialize_metagenome_1(struct _training *);
void initialize_metagenome_2(struct _training *);
void initialize_metagenome_3(struct _training *);
void initialize_metagenome_4(struct _training *);
void initialize_metagenome_5(struct _training *);
void initialize_metagenome_6(struct _training *);
void initialize_metagenome_7(struct _training *);
void initialize_metagenome_8(struct _training *);
void initialize_metagenome_9(struct _training *);
void initialize_metagenome_10(struct _training *);
void initialize_metagenome_11(struct _training *);
void initialize_metagenome_12(struct _training *);
void initialize_metagenome_13(struct _training *);
void initialize_metagenome_14(struct _training *);
void initialize_metagenome_15(struct _training *);
void initialize_metagenome_16(struct _training *);
void initialize_metagenome_17(struct _training *);
void initialize_metagenome_18(struct _training *);
void initialize_metagenome_19(struct _training *);
void initialize_metagenome_20(struct _training *);
void initialize_metagenome_21(struct _training *);
void initialize_metagenome_22(struct _training *);
void initialize_metagenome_23(struct _training *);
void initialize_metagenome_24(struct _training *);
void initialize_metagenome_25(struct _training *);
void initialize_metagenome_26(struct _training *);
void initialize_metagenome_27(struct _training *);
void initialize_metagenome_28(struct _training *);
void initialize_metagenome_29(struct _training *);
void initialize_metagenome_30(struct _training *);
void initialize_metagenome_31(struct _training *);
void initialize_metagenome_32(struct _training *);
void initialize_metagenome_33(struct _training *);
void initialize_metagenome_34(struct _training *);
void initialize_metagenome_35(struct _training *);
void initialize_metagenome_36(struct _training *);
void initialize_metagenome_37(struct _training *);
void initialize_metagenome_38(struct _training *);
void initialize_metagenome_39(struct _training *);
void initialize_metagenome_40(struct _training *);
void initialize_metagenome_41(struct _training *);
void initialize_metagenome_42(struct _training *);
void initialize_metagenome_43(struct _training *);
void initialize_metagenome_44(struct _training *);
void initialize_metagenome_45(struct _training *);
void initialize_metagenome_46(struct _training *);
void initialize_metagenome_47(struct _training *);
void initialize_metagenome_48(struct _training *);
void initialize_metagenome_49(struct _training *);

#endif
