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
  Complete and partial genes are stored separately.  Also records average
  contig length.
******************************************************************************/
void calc_training_set_stats(struct _node *nodes, int initial_node,
                             struct _summary *genome_data, int num_contig,
                             int seq_length)
{
  int path = initial_node;
  int left = 0;
  int right = 0;
  int partial = 0;
  int ctr = 0;

  genome_data->num_complete_genes = 0;
  genome_data->num_partial_genes = 0;
  genome_data->avg_comp_gene_len = 0;
  genome_data->avg_part_gene_len = 0;

  while (path != -1)
  {
    if (nodes[path].status == 0)
    {
      path = nodes[path].trace_forward;
      continue;
    }
    if (nodes[path].strand == 1 && nodes[path].type == START)
    {
      if (nodes[path].cscore > 0.0)
      {
        left = nodes[path].index+1;
      }
      else
      {
        left = -1;
      }
      if (nodes[path].edge == 0)
      {
        partial = 0;
      }
      else
      {
        partial = 1;
      }
    }
    if (nodes[path].strand == -1 && nodes[path].type == STOP)
    {
      left = nodes[path].index-1;
      if (nodes[path].edge == 0)
      {
        partial = 0;
      }
      else
      {
        partial = 1;
      }
    }
    if (nodes[path].strand == 1 && nodes[path].type == STOP)
    {
      right = nodes[path].index+3;
      if (left != -1 && partial == 0 && nodes[path].edge == 0)
      {
        ctr = genome_data->num_complete_genes;
        genome_data->avg_comp_gene_len *= ctr;
        genome_data->avg_comp_gene_len += (right-left+1);
        genome_data->avg_comp_gene_len /= (++ctr);
        genome_data->num_complete_genes++;
      }
      else if (left != -1)
      {
        ctr = genome_data->num_partial_genes;
        genome_data->avg_part_gene_len *= ctr;
        genome_data->avg_part_gene_len += (right-left+1);
        genome_data->avg_part_gene_len /= (++ctr);
        genome_data->num_partial_genes++;
      }
    }
    if (nodes[path].strand == -1 && nodes[path].type == START)
    {
      right = nodes[path].index+1;
      if (nodes[path].cscore > 0.0 && partial == 0 && nodes[path].edge == 0)
      {
        ctr = genome_data->num_complete_genes;
        genome_data->avg_comp_gene_len *= ctr;
        genome_data->avg_comp_gene_len += (right-left+1);
        genome_data->avg_comp_gene_len /= (++ctr);
        genome_data->num_complete_genes++;
      }
      else if (nodes[path].cscore > 0.0)
      {
        ctr = genome_data->num_partial_genes;
        genome_data->avg_part_gene_len *= ctr;
        genome_data->avg_part_gene_len += (right-left+1);
        genome_data->avg_part_gene_len /= (++ctr);
        genome_data->num_partial_genes++;
      }
    }
    path = nodes[path].trace_forward;
  }

  /* Average Contig Length */
  genome_data->num_contig = num_contig;
  genome_data->avg_contig_len = (double)(seq_length - (num_contig-1)*8);
  genome_data->avg_contig_len /= (double)num_contig;
}
