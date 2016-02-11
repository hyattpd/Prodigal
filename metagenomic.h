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

#ifndef SAMPLE_H_
#define SAMPLE_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sequence.h"
#include "training.h"
#include "node.h"

#define NUM_BIN 6
#define NUM_META 50
#define SAMPLE_LEN 120
#define MAX_SAMPLE 200

struct _metagenomic_bin {
  int index;                    /* Index used for sorting */
  int clusnum;                  /* Cluster number */
  char desc[500];               /* Text description of this bin */
  double weight;                /* Current weight/score of this bin */
  double gc;                    /* GC distance from target sequence */
  struct _training *tinf;       /* Pointer to the training file for this bin */
};

void initialize_metagenomic_bins(struct _metagenomic_bin *);
double score_edges(unsigned char *, unsigned char *, int, 
                   struct _training *tinf);
double score_sample(unsigned char *, unsigned char *, int, int, int, struct
                    _training *tinf);
void determine_top_bins(unsigned char *, unsigned char *, int, double,
                        struct _metagenomic_bin *);
int compare_meta_bins(const void *, const void *);

#endif
