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

#ifndef _SUMMARY_H
#define _SUMMARY_H

#include <stdio.h>
#include "node.h"

#define MIN_AVG_TRAIN_GENE_LEN 600
#define MIN_AVG_TRAIN_CTG_LEN 3000

struct _summary {
  int num_contig;
  double avg_contig_len;
  int contig_len_bins[10];
  double avg_contig_gc;
  int contig_gc_bins[12]; /* <25, 25-30, etc., >75 */

  int num_complete_genes;
  double avg_comp_gene_len;
  int comp_gene_len_bins[10];  
  double avg_comp_gene_gc;
  int comp_gene_gc_bins[12];
  int comp_start_bins[5];
  int comp_stop_bins[4];
  int comp_rbs_nosd_bins[4][4][4096];

  int num_partial_genes;
  double avg_part_gene_len;
  int part_gene_len_bins[10];
  double avg_part_gene_gc;
  int part_gene_gc_bins[12];
  int part_start_bins[5];
  int part_stop_bins[4];
  int part_rbs_nosd_bins[4][4][4096];
};

void calc_avg_training_gene_lens(struct _node *, int, struct _summary *);
int bad_train_gene_length(struct _summary);
void calc_avg_train_contig_len(struct _summary *, int);

#endif
