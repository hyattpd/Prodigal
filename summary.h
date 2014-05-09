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
#include "gene.h"

#define MIN_AVG_TRAIN_GENE_LEN 600
#define MIN_AVG_TRAIN_CTG_LEN 3000

void zero_statistics(struct _summary *);
int bad_train_gene_length(struct _summary);
void low_gene_len_warning(int, struct _summary);
void calc_avg_train_contig_len(struct _summary *, int);
void calc_avg_training_gene_lens(struct _node *, int, struct _summary *);

#endif
