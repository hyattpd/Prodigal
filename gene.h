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

#ifndef _GENE_H
#define _GENE_H

#include <stdio.h>
#include "training.h"
#include "node.h"
#include "dprog.h"

#define MAX_GENES 30000

struct _gene {
  int begin;               /* Left end of the gene */
  int end;                 /* Right end of the gene */
  int start_ndx;           /* Index to the start node in the nodes file */
  int stop_ndx;            /* Index to the stop node in the nodes file */
  char gene_data[500];     /* String containing gene information */
  char score_data[500];    /* String containing scoring information */
};

int add_genes(struct _gene *, struct _node *, int);
void record_gene_data(struct _gene *, int, struct _node *, struct _training *,
                      int);
void tweak_final_starts(struct _gene *, int, struct _node *, int, struct
                       _training *);

void print_genes(FILE *, struct _gene *, int, struct _node *, int, int, int,
                 int, char *, struct _training *, char *, char *, char *);
void write_translations(FILE *, struct _gene *, int, struct _node *, 
                        unsigned char *, unsigned char *, unsigned char *, int,
                        struct _training *, int, char *);
void write_nucleotide_seqs(FILE *, struct _gene *, int, struct _node *, 
                           unsigned char *, unsigned char *, unsigned char *,
                           int, struct _training *, int, char *);
double calculate_confidence(double, double);

#endif
