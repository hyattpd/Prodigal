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

#include "summary.h"

/* Memset nodes to 0 and return 0 */
void zero_statistics(struct _summary *statistics)
{
  memset(statistics, 0, sizeof(struct _summary));
}

/******************************************************************************
  Routine to record the average gene length of our initial dynamic programming.
  Complete and partial genes are stored separately.
******************************************************************************/
void calc_avg_training_gene_lens(struct _node *nod, int dbeg, struct _summary
                                 *gstat)
{
  int path, left = 0, right = 0, partial = 0, ctr = 0;

  gstat->num_complete_genes = 0;
  gstat->num_partial_genes = 0;
  gstat->avg_comp_gene_len = 0;
  gstat->avg_part_gene_len = 0;

  if (dbeg == -1)
  {
    return;
  }
  path = dbeg;
  while (nod[path].trace_back != -1)
  {
    path = nod[path].trace_back;
  }
  while (path != -1)
  {
    if (nod[path].eliminate == 1)
    {
      path = nod[path].trace_forward;
      continue;
    }
    if (nod[path].strand == 1 && is_start_node(&nod[path]) == 1)
    {
      if (nod[path].cscore > 0.0)
      {
        left = nod[path].index+1;
      }
      else
      {
        left = -1;
      }
      if (nod[path].edge == 0)
      {
        partial = 0;
      }
      else
      {
        partial = 1;
      }
    }
    if (nod[path].strand == -1 && is_stop_node(&nod[path]) == 1)
    {
      left = nod[path].index-1;
      if (nod[path].edge == 0)
      {
        partial = 0;
      }
      else
      {
        partial = 1;
      }
    }
    if (nod[path].strand == 1 && is_stop_node(&nod[path]) == 1)
    {
      right = nod[path].index+3;
      if (left != -1 && partial == 0 && nod[path].edge == 0)
      {
        ctr = gstat->num_complete_genes;
        gstat->avg_comp_gene_len *= ctr;
        gstat->avg_comp_gene_len += (right-left+1);
        gstat->avg_comp_gene_len /= (++ctr);
        gstat->num_complete_genes++;
      }
      else if (left != -1)
      {
        ctr = gstat->num_partial_genes;
        gstat->avg_part_gene_len *= ctr;
        gstat->avg_part_gene_len += (right-left+1);
        gstat->avg_part_gene_len /= (++ctr);
        gstat->num_partial_genes++;
      }
    }
    if (nod[path].strand == -1 && is_start_node(&nod[path]) == 1)
    {
      right = nod[path].index+1;
      if (nod[path].cscore > 0.0 && partial == 0 && nod[path].edge == 0)
      {
        ctr = gstat->num_complete_genes;
        gstat->avg_comp_gene_len *= ctr;
        gstat->avg_comp_gene_len += (right-left+1);
        gstat->avg_comp_gene_len /= (++ctr);
        gstat->num_complete_genes++;
      }
      else if (nod[path].cscore > 0.0)
      {
        ctr = gstat->num_partial_genes;
        gstat->avg_part_gene_len *= ctr;
        gstat->avg_part_gene_len += (right-left+1);
        gstat->avg_part_gene_len /= (++ctr);
        gstat->num_partial_genes++;
      }
    }
    path = nod[path].trace_forward;
  }
}

/******************************************************************************
  Look at average complete gene length.  If it's above the minimum acceptable
  length, return 0.  If it's too small, see if this is due to the sequence
  being in tons of contigs.  If so, return 1.  If we still have no explanation
  for why it's low, return 2.
******************************************************************************/
int bad_train_gene_length(struct _summary gstat)
{
  if (gstat.avg_comp_gene_len > MIN_AVG_TRAIN_GENE_LEN)
  {
    return 0;
  }
  if (gstat.avg_contig_len < MIN_AVG_TRAIN_CTG_LEN ||
      gstat.num_partial_genes > gstat.num_complete_genes)
  {
    return 1;
  }
  return 2;
}

/* Output a warning for low average gene length */
void low_gene_len_warning(int flag, struct _summary gstat)
{
  if (flag < 2)
  {
    fprintf(stderr, "\nWarning: training sequence is highly fragmented.\n");
    fprintf(stderr, "You may get better results with the ");
    fprintf(stderr, "'-m anon' option.\n\n");
  }
  else
  {
    fprintf(stderr, "\nWarning: Average training gene length is");
    fprintf(stderr, " low (%.1f).\n", gstat.avg_comp_gene_len);
    fprintf(stderr, "Double check translation table or check for");
    fprintf(stderr, " pseudogenes/gene decay.\n\n");
  }
}

/******************************************************************************
  Calculate average length of contigs in training set.  Since we added gaps in
  between each contig, we subtract 8 bases per contig to get the true length.
******************************************************************************/
void calc_avg_train_contig_len(struct _summary *gstat, int slen)
{
  gstat->avg_contig_len = (double)(slen - (gstat->num_contig-1)*8);
  gstat->avg_contig_len /= (double)gstat->num_contig;
}
