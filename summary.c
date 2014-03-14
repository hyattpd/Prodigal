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

#include "summary.h"

/*******************************************************************************
  Routine to record the average gene length of our initial dynamic programming.
  Complete and partial genes are stored separately.
*******************************************************************************/
void calc_avg_training_gene_lens(struct _node *nod, int dbeg, struct _summary
                                 *gstat) {
  int path, left = 0, right = 0, partial = 0, ctr = 0;

  if(dbeg == -1) return;
  path = dbeg;
  while(nod[path].traceb != -1) path = nod[path].traceb;
  while(path != -1) {
    if(nod[path].elim == 1) { path = nod[path].tracef; continue; }
    if(nod[path].strand == 1 && nod[path].type != STOP) {
      if(nod[path].cscore > 0.0) left = nod[path].ndx+1;
      else left = -1;
      if(nod[path].edge == 0) partial = 0;
      else partial = 1;
    }
    if(nod[path].strand == -1 && nod[path].type == STOP) {
      left = nod[path].ndx-1;
      if(nod[path].edge == 0) partial = 0;
      else partial = 1;
    }
    if(nod[path].strand == 1 && nod[path].type == STOP) {
      right = nod[path].ndx+3;
      if(left != -1 && partial == 0 && nod[path].edge == 0) {
        ctr = gstat->num_complete_genes;
        gstat->avg_comp_gene_len *= ctr;
        gstat->avg_comp_gene_len += (right-left+1);
        gstat->avg_comp_gene_len /= (++ctr);
        gstat->num_complete_genes++;
      }
      else if(left != -1) {
        ctr = gstat->num_partial_genes;
        gstat->avg_part_gene_len *= ctr;
        gstat->avg_part_gene_len += (right-left+1);
        gstat->avg_part_gene_len /= (++ctr);
        gstat->num_partial_genes++;
      }
    }
    if(nod[path].strand == -1 && nod[path].type != STOP) {
      right = nod[path].ndx+1;
      if(nod[path].cscore > 0.0 && partial == 0 && nod[path].edge == 0) {
        ctr = gstat->num_complete_genes;
        gstat->avg_comp_gene_len *= ctr;
        gstat->avg_comp_gene_len += (right-left+1);
        gstat->avg_comp_gene_len /= (++ctr);
        gstat->num_complete_genes++;
      }
      else if(nod[path].cscore > 0.0) {
        ctr = gstat->num_partial_genes;
        gstat->avg_part_gene_len *= ctr;
        gstat->avg_part_gene_len += (right-left+1);
        gstat->avg_part_gene_len /= (++ctr);
        gstat->num_partial_genes++;
      }
    }
    path = nod[path].tracef;
  }
}

/******************************************************************************
  Look at average complete gene length.  If it's above the minimum acceptable
  length, return 0.  If it's too small, see if this is due to the sequence 
  being in tons of contigs.  If so, return 1.  If we still have no explanation
  for why it's low, return 2.
******************************************************************************/
int bad_train_gene_length(struct _summary gstat) {
  if(gstat.avg_comp_gene_len > MIN_AVG_TRAIN_GENE_LEN) return 0;
  if(gstat.avg_contig_len < MIN_AVG_TRAIN_CTG_LEN ||
     gstat.num_partial_genes > gstat.num_complete_genes) return 1;
  return 2;
}

void calc_avg_train_contig_len(struct _summary *gstat, int slen) {
  gstat->avg_contig_len = (double)slen/(double)gstat->num_contig;
}
