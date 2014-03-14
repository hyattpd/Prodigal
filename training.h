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

#ifndef _TRAIN_H
#define _TRAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datatypes.h"
#include "node.h"

int write_training_file(char *, struct _training *);
int read_training_file(char *, struct _training *);

void record_gc_frame_bias(struct _training *, unsigned char *, int, struct 
                          _node *, int);
void calc_dicodon_gene(struct _training *, unsigned char *, unsigned char *,
                       int, struct _node *, int);

void train_starts_sd(unsigned char *, unsigned char *, int, struct _node *,
                        int, struct _training *);
void train_starts_nonsd(unsigned char *, unsigned char *, int, struct _node *,
                        int, struct _training *);

void count_upstream_composition(unsigned char *, int, int, int,
                                struct _training *);
void build_coverage_map(double [4][4][4096], int [4][4][4096], double, int);
void update_motif_counts(double [4][4][4096], double *, unsigned char *,
                         unsigned char *, int, struct _node *, int);

void determine_sd_usage(struct _training *);

void initialize_preset_genome_0(struct _training *);
void initialize_preset_genome_1(struct _training *);
void initialize_preset_genome_2(struct _training *);
void initialize_preset_genome_3(struct _training *);
void initialize_preset_genome_4(struct _training *);
void initialize_preset_genome_5(struct _training *);
void initialize_preset_genome_6(struct _training *);
void initialize_preset_genome_7(struct _training *);
void initialize_preset_genome_8(struct _training *);
void initialize_preset_genome_9(struct _training *);
void initialize_preset_genome_10(struct _training *);
void initialize_preset_genome_11(struct _training *);
void initialize_preset_genome_12(struct _training *);
void initialize_preset_genome_13(struct _training *);
void initialize_preset_genome_14(struct _training *);
void initialize_preset_genome_15(struct _training *);
void initialize_preset_genome_16(struct _training *);
void initialize_preset_genome_17(struct _training *);
void initialize_preset_genome_18(struct _training *);
void initialize_preset_genome_19(struct _training *);
void initialize_preset_genome_20(struct _training *);
void initialize_preset_genome_21(struct _training *);
void initialize_preset_genome_22(struct _training *);
void initialize_preset_genome_23(struct _training *);
void initialize_preset_genome_24(struct _training *);
void initialize_preset_genome_25(struct _training *);
void initialize_preset_genome_26(struct _training *);
void initialize_preset_genome_27(struct _training *);
void initialize_preset_genome_28(struct _training *);
void initialize_preset_genome_29(struct _training *);
void initialize_preset_genome_30(struct _training *);
void initialize_preset_genome_31(struct _training *);
void initialize_preset_genome_32(struct _training *);
void initialize_preset_genome_33(struct _training *);
void initialize_preset_genome_34(struct _training *);
void initialize_preset_genome_35(struct _training *);
void initialize_preset_genome_36(struct _training *);
void initialize_preset_genome_37(struct _training *);
void initialize_preset_genome_38(struct _training *);
void initialize_preset_genome_39(struct _training *);
void initialize_preset_genome_40(struct _training *);
void initialize_preset_genome_41(struct _training *);
void initialize_preset_genome_42(struct _training *);
void initialize_preset_genome_43(struct _training *);
void initialize_preset_genome_44(struct _training *);
void initialize_preset_genome_45(struct _training *);
void initialize_preset_genome_46(struct _training *);
void initialize_preset_genome_47(struct _training *);
void initialize_preset_genome_48(struct _training *);
void initialize_preset_genome_49(struct _training *);

#endif
