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

#include "training.h"

/* Reads a training file to use for gene prediction */
int read_training_file(char *fn, struct _training *tinf)
{
  size_t rv;
  FILE *fh;
  fh = fopen(fn, "rb");
  if (fh == NULL)
  {
    return -1;
  }
  rv = fread(tinf, sizeof(struct _training), 1, fh);
  fclose(fh);
  if (rv != 1)
  {
    return -1;
  }
  return 0;
}

/* Writes a training file to use for a later run of gene prediction */
int write_training_file(char *fn, struct _training *tinf)
{
  size_t rv;
  FILE *fh;
  fh = fopen(fn, "wb");
  if (fh == NULL)
  {
    return -1;
  }
  rv = fwrite(tinf, sizeof(struct _training), 1, fh);
  fclose(fh);
  if (rv != 1)
  {
    return -1;
  }
  return 0;
}

/* Build the training set using the supplied genetic code */
void build_training_set(struct _node *nodes, struct _training *tinf, struct
                        _summary *statistics, unsigned char *seq, unsigned
                        char *rseq, unsigned char *useq, int slen, int
                        *num_node, int closed, int cross_gaps)
{
  int nn = 0, ipath = 0;

  /***********************************************************************
    Find all the potential starts and stops, sort them, and create a
    comprehensive list of nodes for dynamic programming.
  ***********************************************************************/
  if (*num_node > 0)
  {
    zero_nodes(nodes, *num_node);
  }
  nn = add_nodes(seq, rseq, useq, slen, nodes, closed, cross_gaps,
                 tinf->trans_table);
  *num_node = nn;
  qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

  /***********************************************************************
    Scan all the ORFS looking for a potential GC bias in a particular
    codon position.  This information will be used to acquire a good
    initial set of genes.
  ***********************************************************************/
  record_gc_frame_bias(tinf, seq, slen, nodes, nn);

  /***********************************************************************
    Do an initial dynamic programming routine with just the GC frame
    bias used as a scoring function.  This will get an initial set of
    genes to train on.
  ***********************************************************************/
  record_overlapping_starts(nodes, nn, tinf->start_weight, 0);
  ipath = dynamic_programming(nodes, nn, tinf->bias, tinf->start_weight, 0);

  /***********************************************************************
    Gather dicodon statistics for the training set.  Score the entire set
    of nodes.
  ***********************************************************************/
  calc_dicodon_gene(tinf, seq, rseq, slen, nodes, ipath);
  raw_coding_score(seq, rseq, slen, nodes, nn, tinf->trans_table, tinf->gc,
                   tinf->gene_dc);

  /***********************************************************************
    Gather statistics about average gene length to see if the training
    set looks good.
  ***********************************************************************/
  calc_avg_training_gene_lens(nodes, ipath, statistics);
}

/* Records the GC frame bias from the node GC statistics */
void record_gc_frame_bias(struct _training *tinf, unsigned char *seq, int slen,
                          struct _node *nod, int nn)
{
  int i, len, *gc_frame;
  double tot;

  gc_frame = calc_most_gc_frame(seq, slen);
  frame_score(gc_frame, nod, nn);
  free(gc_frame);

  for (i = 0; i < 3; i++)
  {
    tinf->bias[i] = 0.0;
  }
  for (i = 0; i < nn; i++)
  {
    if (is_start_node(&nod[i]) == 1)
    {
      len = abs(nod[i].stop_val-nod[i].index)+1;
      tinf->bias[nod[i].gc_bias] +=
        (nod[i].gc_score[nod[i].gc_bias]*len)/1000.0;
    }
  }
  tot = tinf->bias[0] + tinf->bias[1] + tinf->bias[2];
  for (i = 0; i < 3; i++)
  {
    tinf->bias[i] *= (3.0/tot);
  }
}

/******************************************************************************
  Simple routine that calculates the dicodon frequency in genes and in the
  background, and then stores the log likelihood of each 6-mer relative to the
  background.
******************************************************************************/

void calc_dicodon_gene(struct _training *tinf, unsigned char *seq, unsigned
                       char *rseq, int slen, struct _node *nod, int dbeg)
{
  int i, path, counts[4096], glob = 0;
  int left, right, in_gene;
  double prob[4096], bg[4096];

  for (i = 0; i < 4096; i++)
  {
    counts[i] = 0;
    prob[i] = 0.0;
    bg[i] = 0.0;
  }
  left = -1;
  right = -1;
  calc_mer_bg(6, seq, rseq, slen, bg);
  path = dbeg;
  in_gene = 0;
  while (path != -1)
  {
    if (nod[path].strand == -1 && is_start_node(&nod[path]) == 1)
    {
      in_gene = -1;
      left = slen-nod[path].index-1;
    }
    if (nod[path].strand == 1 && is_stop_node(&nod[path]) == 1)
    {
      in_gene = 1;
      right = nod[path].index+2;
    }
    if (in_gene == -1 && nod[path].strand == -1 &&
        is_stop_node(&nod[path]) == 1)
    {
      right = slen-nod[path].index+1;
      for (i = left; i < right-5; i+=3)
      {
        counts[mer_index(6, rseq, i)]++;
        glob++;
      }
      in_gene = 0;
    }
    if (in_gene == 1 && nod[path].strand == 1 &&
        is_start_node(&nod[path]) == 1)
    {
      left = nod[path].index;
      for (i = left; i < right-5; i+=3)
      {
        counts[mer_index(6, seq, i)]++;
        glob++;
      }
      in_gene = 0;
    }
    path = nod[path].trace_back;
  }
  for (i = 0; i < 4096; i++)
  {
    prob[i] = (counts[i]*1.0)/(glob*1.0);
    if (prob[i] == 0 && bg[i] != 0)
    {
      tinf->gene_dc[i] = -5.0;
    }
    else if (bg[i] == 0)
    {
      tinf->gene_dc[i] = 0.0;
    }
    else
    {
      tinf->gene_dc[i] = log(prob[i]/bg[i]);
    }
    if (tinf->gene_dc[i] > 5.0)
    {
      tinf->gene_dc[i] = 5.0;
    }
    if (tinf->gene_dc[i] < -5.0)
    {
      tinf->gene_dc[i] = -5.0;
    }
  }
}

/******************************************************************************
  Iterative Algorithm to train starts.  It begins with all the highest coding
  starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts moving
  starts around attempting to match these discoveries.  This start trainer is
  for Shine-Dalgarno motifs only.
******************************************************************************/
void train_starts_sd(unsigned char *seq, unsigned char *rseq, int slen,
                  struct _node *nod, int nn, struct _training *tinf)
{
  int i, j, fr, rbs[3], type[3], bindex[3], max_rb;
  double sum, wt, rbg[28], rreal[28], best[3], sthresh = 35.0;
  double tbg[3], treal[3];

  wt = tinf->start_weight;
  for (j = 0; j < 3; j++)
  {
    tinf->type_wt[j] = 0.0;
  }
  for (j = 0; j < 28; j++)
  {
    tinf->rbs_wt[j] = 0.0;
  }
  for (i = 0; i < 32; i++)
  {
    for (j = 0; j < 4; j++)
    {
      tinf->ups_comp[i][j] = 0.0;
    }
  }

  /* Build the background of random types */
  for (i = 0; i < 3; i++)
  {
    tbg[i] = 0.0;
  }
  for (i = 0; i < nn; i++)
  {
    if (is_stop_node(&nod[i]) == 1)
    {
      continue;
    }
    tbg[nod[i].type] += 1.0;
  }
  sum = 0.0;
  for (i = 0; i < 3; i++)
  {
    sum += tbg[i];
  }
  for (i = 0; i < 3; i++)
  {
    tbg[i] /= sum;
  }

  /* Iterate 10 times through the list of nodes                         */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few     */
  /* extra to be safe)                                                  */
  for (i = 0; i < 10; i++)
  {

    /* Recalculate the RBS motif background */
    for (j = 0; j < 28; j++)
    {
      rbg[j] = 0.0;
    }
    for (j = 0; j < nn; j++)
    {
      if (is_stop_node(&nod[j]) == 1 || nod[j].edge == 1)
      {
        continue;
      }
      if (tinf->rbs_wt[nod[j].rbs[0]] > tinf->rbs_wt[nod[j].rbs[1]]+1.0 ||
          nod[j].rbs[1] == 0)
      {
        max_rb = nod[j].rbs[0];
      }
      else if (tinf->rbs_wt[nod[j].rbs[0]] < tinf->rbs_wt[nod[j].rbs[1]]-1.0 ||
               nod[j].rbs[0] == 0)
      {
        max_rb = nod[j].rbs[1];
      }
      else
      {
        max_rb = (int)dmax(nod[j].rbs[0], nod[j].rbs[1]);
      }
      rbg[max_rb] += 1.0;
    }
    sum = 0.0;
    for (j = 0; j < 28; j++)
    {
      sum += rbg[j];
    }
    for (j = 0; j < 28; j++)
    {
      rbg[j] /= sum;
    }

    for (j = 0; j < 28; j++)
    {
      rreal[j] = 0.0;
    }
    for (j = 0; j < 3; j++)
    {
      treal[j] = 0.0;
    }

    /* Forward strand pass */
    for (j = 0; j < 3; j++)
    {
      best[j] = 0.0;
      bindex[j] = -1;
      rbs[j] = 0;
      type[j] = 0;
    }
    for (j = 0; j < nn; j++)
    {
      if (is_start_node(&nod[j]) == 1 && nod[j].edge == 1)
      {
        continue;
      }
      fr = (nod[j].index)%3;
      if (is_stop_node(&nod[j]) == 1 && nod[j].strand == 1)
      {
        if (best[fr] >= sthresh && nod[bindex[fr]].index%3 == fr)
        {
          rreal[rbs[fr]] += 1.0;
          treal[type[fr]] += 1.0;
          if (i == 9)
          {
            count_upstream_composition(seq, slen, 1, nod[bindex[fr]].index,
                                       tinf);
          }
        }
        best[fr] = 0.0;
        bindex[fr] = -1;
        rbs[fr] = 0;
        type[fr] = 0;
      }
      else if (nod[j].strand == 1)
      {
        if (tinf->rbs_wt[nod[j].rbs[0]] > tinf->rbs_wt[nod[j].rbs[1]]+1.0 ||
            nod[j].rbs[1] == 0)
        {
          max_rb = nod[j].rbs[0];
        }
        else if (tinf->rbs_wt[nod[j].rbs[0]] < tinf->rbs_wt[nod[j].rbs[1]]-1.0
                 || nod[j].rbs[0] == 0)
        {
          max_rb = nod[j].rbs[1];
        }
        else
        {
          max_rb = (int)dmax(nod[j].rbs[0], nod[j].rbs[1]);
        }
        if (nod[j].cscore + wt*tinf->rbs_wt[max_rb] +
            wt*tinf->type_wt[nod[j].type] >= best[fr])
        {
          best[fr] = nod[j].cscore + wt*tinf->rbs_wt[max_rb];
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bindex[fr] = j;
          type[fr] = nod[j].type;
          rbs[fr] = max_rb;
        }
      }
    }

    /* Reverse strand pass */
    for (j = 0; j < 3; j++)
    {
      best[j] = 0.0;
      bindex[j] = -1;
      rbs[j] = 0;
      type[j] = 0;
    }
    for (j = nn-1; j >= 0; j--)
    {
      if (is_start_node(&nod[j]) == 1 && nod[j].edge == 1)
      {
        continue;
      }
      fr = (nod[j].index)%3;
      if (is_stop_node(&nod[j]) == 1 && nod[j].strand == -1)
      {
        if (best[fr] >= sthresh && nod[bindex[fr]].index%3 == fr)
        {
          rreal[rbs[fr]] += 1.0;
          treal[type[fr]] += 1.0;
          if (i == 9)
          {
            count_upstream_composition(rseq, slen, -1, nod[bindex[fr]].index,
                                       tinf);
          }
        }
        best[fr] = 0.0;
        bindex[fr] = -1;
        rbs[fr] = 0;
        type[fr] = 0;
      }
      else if (nod[j].strand == -1)
      {
        if (tinf->rbs_wt[nod[j].rbs[0]] > tinf->rbs_wt[nod[j].rbs[1]]+1.0 ||
            nod[j].rbs[1] == 0)
        {
          max_rb = nod[j].rbs[0];
        }
        else if (tinf->rbs_wt[nod[j].rbs[0]] < tinf->rbs_wt[nod[j].rbs[1]]-1.0
                 || nod[j].rbs[0] == 0)
        {
          max_rb = nod[j].rbs[1];
        }
        else
        {
          max_rb = (int)dmax(nod[j].rbs[0], nod[j].rbs[1]);
        }
        if (nod[j].cscore + wt*tinf->rbs_wt[max_rb] +
            wt*tinf->type_wt[nod[j].type] >= best[fr])
        {
          best[fr] = nod[j].cscore + wt*tinf->rbs_wt[max_rb];
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bindex[fr] = j;
          type[fr] = nod[j].type;
          rbs[fr] = max_rb;
        }
      }
    }

    sum = 0.0;
    for (j = 0; j < 28; j++)
    {
      sum += rreal[j];
    }
    if (sum == 0.0)
    {
      for (j = 0; j < 28; j++)
      {
        tinf->rbs_wt[j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 28; j++)
      {
        rreal[j] /= sum;
        if (rbg[j] != 0)
        {
          tinf->rbs_wt[j] = log(rreal[j]/rbg[j]);
        }
        else
        {
          tinf->rbs_wt[j] = -4.0;
        }
        if (tinf->rbs_wt[j] > 4.0)
        {
          tinf->rbs_wt[j] = 4.0;
        }
        if (tinf->rbs_wt[j] < -4.0)
        {
          tinf->rbs_wt[j] = -4.0;
        }
      }
    }
    sum = 0.0;
    for (j = 0; j < 3; j++)
    {
      sum += treal[j];
    }
    if (sum == 0.0)
    {
      for (j = 0; j < 3; j++)
      {
        tinf->type_wt[j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 3; j++)
      {
        treal[j] /= sum;
        if (tbg[j] != 0)
        {
          tinf->type_wt[j] = log(treal[j]/tbg[j]);
        }
        else
        {
          tinf->type_wt[j] = -4.0;
        }
        if (tinf->type_wt[j] > 4.0)
        {
          tinf->type_wt[j] = 4.0;
        }
        if (tinf->type_wt[j] < -4.0)
        {
          tinf->type_wt[j] = -4.0;
        }
      }
    }
    if (sum <= (double)nn/2000.0)
    {
      sthresh /= 2.0;
    }
  }

  /* Convert upstream base composition to a log score */
  for (i = 0; i < 32; i++)
  {
    sum = 0.0;
    for (j = 0; j < 4; j++)
    {
      sum += tinf->ups_comp[i][j];
    }
    if (sum == 0.0)
    {
      for (j = 0; j < 4; j++)
      {
        tinf->ups_comp[i][j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 4; j++)
      {
        tinf->ups_comp[i][j] /= sum;
        if (tinf->gc > 0.1 && tinf->gc < 0.9)
        {
          if (j == 0 || j == 3)
          {
            tinf->ups_comp[i][j] =
              log(tinf->ups_comp[i][j]*2.0/(1.0-tinf->gc));
          }
          else
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/tinf->gc);
          }
        }
        else if (tinf->gc <= 0.1)
        {
          if (j == 0 || j == 3)
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
          }
          else
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
          }
        }
        else
        {
          if (j == 0 || j == 3)
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
          }
          else
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
          }
        }
        if (tinf->ups_comp[i][j] > 4.0)
        {
          tinf->ups_comp[i][j] = 4.0;
        }
        if (tinf->ups_comp[i][j] < -4.0)
        {
          tinf->ups_comp[i][j] = -4.0;
        }
      }
    }
  }
}

/******************************************************************************
  Iterative Algorithm to train starts.  It begins with all the highest coding
  starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts moving
  starts around attempting to match these discoveries.  Unlike the SD
  algorithm, it allows for any popular motif to be discovered.
******************************************************************************/
void train_starts_nonsd(unsigned char *seq, unsigned char *rseq, int slen,
                  struct _node *nod, int nn, struct _training *tinf)
{
  int i, j, k, l, fr, bindex[3], mgood[4][4][4096], stage;
  double sum, ngenes, wt = tinf->start_weight, best[3], sthresh = 35.0;
  double tbg[3], treal[3];
  double mbg[4][4][4096], mreal[4][4][4096], zbg, zreal;

  for (i = 0; i < 32; i++)
  {
    for (j = 0; j < 4; j++)
    {
      tinf->ups_comp[i][j] = 0.0;
    }
  }

  /* Build the background of random types */
  for (i = 0; i < 3; i++)
  {
    tinf->type_wt[i] = 0.0;
  }
  for (i = 0; i < 3; i++)
  {
    tbg[i] = 0.0;
  }
  for (i = 0; i < nn; i++)
  {
    if (is_stop_node(&nod[i]) == 1)
    {
      continue;
    }
    tbg[nod[i].type] += 1.0;
  }
  sum = 0.0;
  for (i = 0; i < 3; i++)
  {
    sum += tbg[i];
  }
  for (i = 0; i < 3; i++)
  {
    tbg[i] /= sum;
  }

  /* Iterate 20 times through the list of nodes                         */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few      */
  /* extra to be safe)                                                  */
  for (i = 0; i < 20; i++)
  {

    /* Determine which stage of motif finding we're in */
    if (i < 4)
    {
      stage = 0;
    }
    else if (i < 12)
    {
      stage = 1;
    }
    else
    {
      stage = 2;
    }

    /* Recalculate the upstream motif background and set 'real' counts to 0 */
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          mbg[j][k][l] = 0.0;
        }
      }
    }
    zbg = 0.0;
    for (j = 0; j < nn; j++)
    {
      if (is_stop_node(&nod[j]) == 1 || nod[j].edge == 1)
      {
        continue;
      }
      find_best_upstream_motif(tinf, seq, rseq, slen, &nod[j], stage);
      update_motif_counts(mbg, &zbg, seq, rseq, slen, &(nod[j]), stage);
    }
    sum = 0.0;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          sum += mbg[j][k][l];
        }
      }
    }
    sum += zbg;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          mbg[j][k][l] /= sum;
        }
      }
    }
    zbg /= sum;

    /* Reset counts of 'real' motifs/types to 0 */
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          mreal[j][k][l] = 0.0;
        }
      }
    }
    zreal = 0.0;
    for (j = 0; j < 3; j++)
    {
      treal[j] = 0.0;
    }
    ngenes = 0.0;

    /* Forward strand pass */
    for (j = 0; j < 3; j++)
    {
      best[j] = 0.0;
      bindex[j] = -1;
    }
    for (j = 0; j < nn; j++)
    {
      if (is_start_node(&nod[j]) == 1 && nod[j].edge == 1)
      {
        continue;
      }
      fr = (nod[j].index)%3;
      if (is_stop_node(&nod[j]) == 1 && nod[j].strand == 1)
      {
        if (best[fr] >= sthresh)
        {
          ngenes += 1.0;
          treal[nod[bindex[fr]].type] += 1.0;
          update_motif_counts(mreal, &zreal, seq, rseq, slen,
                              &(nod[bindex[fr]]), stage);
          if (i == 19)
          {
            count_upstream_composition(seq, slen, 1, nod[bindex[fr]].index,
                                       tinf);
          }
        }
        best[fr] = 0.0;
        bindex[fr] = -1;
      }
      else if (nod[j].strand == 1)
      {
        if (nod[j].cscore + wt*nod[j].mot.score + wt*tinf->type_wt[nod[j].type]
            >= best[fr])
        {
          best[fr] = nod[j].cscore + wt*nod[j].mot.score;
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bindex[fr] = j;
        }
      }
    }

    /* Reverse strand pass */
    for (j = 0; j < 3; j++)
    {
      best[j] = 0.0;
      bindex[j] = -1;
    }
    for (j = nn-1; j >= 0; j--)
    {
      if (is_start_node(&nod[j]) == 1 && nod[j].edge == 1)
      {
        continue;
      }
      fr = (nod[j].index)%3;
      if (is_stop_node(&nod[j]) == 1 && nod[j].strand == -1)
      {
        if (best[fr] >= sthresh)
        {
          ngenes += 1.0;
          treal[nod[bindex[fr]].type] += 1.0;
          update_motif_counts(mreal, &zreal, seq, rseq, slen,
                              &(nod[bindex[fr]]), stage);
          if (i == 19)
          {
            count_upstream_composition(rseq, slen, -1, nod[bindex[fr]].index,
                                       tinf);
          }
        }
        best[fr] = 0.0;
        bindex[fr] = -1;
      }
      else if (nod[j].strand == -1)
      {
        if (nod[j].cscore + wt*nod[j].mot.score + wt*tinf->type_wt[nod[j].type]
            >= best[fr])
        {
          best[fr] = nod[j].cscore + wt*nod[j].mot.score;
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bindex[fr] = j;
        }
      }
    }

    /* Update the log likelihood weights for type and RBS motifs */
    if (stage < 2)
    {
      build_coverage_map(mreal, mgood, ngenes, stage);
    }
    sum = 0.0;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          sum += mreal[j][k][l];
        }
      }
    }
    sum += zreal;
    if (sum == 0.0)
    {
      for (j = 0; j < 4; j++)
      {
        for (k = 0; k < 4; k++)
        {
          for (l = 0; l < 4096; l++)
          {
            tinf->mot_wt[j][k][l] = 0.0;
          }
        }
      }
      tinf->no_mot = 0.0;
    }
    else
    {
      for (j = 0; j < 4; j++)
      {
        for (k = 0; k < 4; k++)
        {
          for (l = 0; l < 4096; l++)
          {
            if (mgood[j][k][l] == 0)
            {
              zreal += mreal[j][k][l];
              zbg += mreal[j][k][l];
              mreal[j][k][l] = 0.0;
              mbg[j][k][l] = 0.0;
            }
            mreal[j][k][l] /= sum;
            if (mbg[j][k][l] != 0)
            {
              tinf->mot_wt[j][k][l] = log(mreal[j][k][l]/mbg[j][k][l]);
            }
            else
            {
              tinf->mot_wt[j][k][l] = -4.0;
            }
            if (tinf->mot_wt[j][k][l] > 4.0)
            {
              tinf->mot_wt[j][k][l] = 4.0;
            }
            if (tinf->mot_wt[j][k][l] < -4.0)
            {
              tinf->mot_wt[j][k][l] = -4.0;
            }
          }
        }
      }
    }
    zreal /= sum;
    if (zbg != 0)
    {
      tinf->no_mot = log(zreal/zbg);
    }
    else
    {
      tinf->no_mot = -4.0;
    }
    if (tinf->no_mot > 4.0)
    {
      tinf->no_mot = 4.0;
    }
    if (tinf->no_mot < -4.0)
    {
      tinf->no_mot = -4.0;
    }
    sum = 0.0;
    for (j = 0; j < 3; j++)
    {
      sum += treal[j];
    }
    if (sum == 0.0)
    {
      for (j = 0; j < 3; j++)
      {
        tinf->type_wt[j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 3; j++)
      {
        treal[j] /= sum;
        if (tbg[j] != 0)
        {
          tinf->type_wt[j] = log(treal[j]/tbg[j]);
        }
        else
        {
          tinf->type_wt[j] = -4.0;
        }
        if (tinf->type_wt[j] > 4.0)
        {
          tinf->type_wt[j] = 4.0;
        }
        if (tinf->type_wt[j] < -4.0)
        {
          tinf->type_wt[j] = -4.0;
        }
      }
    }
    if (sum <= (double)nn/2000.0)
    {
      sthresh /= 2.0;
    }
  }
  /* Convert upstream base composition to a log score */
  for (i = 0; i < 32; i++)
  {
    sum = 0.0;
    for (j = 0; j < 4; j++)
    {
      sum += tinf->ups_comp[i][j];
    }
    if (sum == 0.0)
    {
      for (j = 0; j < 4; j++)
      {
        tinf->ups_comp[i][j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 4; j++)
      {
        tinf->ups_comp[i][j] /= sum;
        if (tinf->gc > 0.1 && tinf->gc < 0.9)
        {
          if (j == 0 || j == 3)
          {
            tinf->ups_comp[i][j] =
              log(tinf->ups_comp[i][j]*2.0/(1.0-tinf->gc));
          }
          else
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/tinf->gc);
          }
        }
        else if (tinf->gc <= 0.1)
        {
          if (j == 0 || j == 3)
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
          }
          else
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
          }
        }
        else
        {
          if (j == 0 || j == 3)
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
          }
          else
          {
            tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
          }
        }
        if (tinf->ups_comp[i][j] > 4.0)
        {
          tinf->ups_comp[i][j] = 4.0;
        }
        if (tinf->ups_comp[i][j] < -4.0)
        {
          tinf->ups_comp[i][j] = -4.0;
        }
      }
    }
  }
}

/******************************************************************************
  For a given start, record the base composition of the upstream region at
  positions -1 and -2 and -15 to -44.  This will be used to supplement the
  SD (or other) motif finder with additional information.
******************************************************************************/
void count_upstream_composition(unsigned char *seq, int slen, int strand,
                                int pos, struct _training *tinf)
{
  int i, start, count = 0;
  if (strand == 1)
  {
    start = pos;
  }
  else
  {
    start = slen-1-pos;
  }

  for (i = 1; i < 45; i++)
  {
    if (i > 2 && i < 15)
    {
      continue;
    }
    tinf->ups_comp[count][mer_index(1, seq, start-i)]++;
    count++;
  }
}

/******************************************************************************
  Update the motif counts from a putative "real" start.  This is done in three
  stages.  In stage 0, all motifs sizes 3-6bp in the region with spacer 3-15bp
  are counted.  In stage 1, only the best motif and all its subsets are
  counted (e.g. for AGGAG, we would count AGGAG, AGGA, GGAG, AGG, GGA, and
  GAG).  In stage 2, only the best single motif is counted.
******************************************************************************/
void update_motif_counts(double mcnt[4][4][4096], double *zero, unsigned char
                         *seq, unsigned char *rseq, int slen,
                         struct _node *nod, int stage)
{
  int i, j, k, start, spaceindex;
  unsigned char *wseq;
  struct _motif *mot = &(nod->mot);

  if (is_stop_node(nod) == 1 || nod->edge == 1)
  {
    return;
  }
  if (mot->len == 0)
  {
    *zero += 1.0;
    return;
  }

  if (nod->strand == 1)
  {
    wseq = seq;
    start = nod->index;
  }
  else
  {
    wseq = rseq;
    start = slen-1-nod->index;
  }

  /* Stage 0:  Count all motifs.  If a motif is detected, */
  /* it is counted for every distance in stage 0.  This   */
  /* is done to make sure off-distance good motifs are    */
  /* recognized.                                          */
  if (stage == 0)
  {
    for (i = 3; i >= 0; i--)
    {
      for (j = start-18-i; j <= start-6-i; j++)
      {
        if (j < 0)
        {
          continue;
        }
        if (j <= start-16-i)
        {
          spaceindex = 3;
        }
        else if (j <= start-14-i)
        {
          spaceindex = 2;
        }
        else if (j >= start-7-i)
        {
          spaceindex = 1;
        }
        else
        {
          spaceindex = 0;
        }
        for (k = 0; k < 4; k++)
        {
          mcnt[i][k][mer_index(i+3, wseq, j)] += 1.0;
        }
      }
    }
  }
  /* Stage 1:  Count only the best motif, but also count  */
  /* all its sub-motifs.                                  */
  else if (stage == 1)
  {
    mcnt[mot->len-3][mot->spaceindex][mot->index] += 1.0;
    for (i = 0; i < mot->len-3; i++)
    {
      for (j = start-(mot->spacer)-(mot->len); j <= start-(mot->spacer)-(i+3);
           j++)
      {
        if (j < 0)
        {
          continue;
        }
        if (j <= start-16-i)
        {
          spaceindex = 3;
        }
        else if (j <= start-14-i)
        {
          spaceindex = 2;
        }
        else if (j >= start-7-i)
        {
          spaceindex = 1;
        }
        else
        {
          spaceindex = 0;
        }
        mcnt[i][spaceindex][mer_index(i+3, wseq, j)] += 1.0;
      }
    }
  }
  /* Stage 2:  Only count the highest scoring motif. */
  else if (stage == 2)
  {
    mcnt[mot->len-3][mot->spaceindex][mot->index] += 1.0;
  }
}

/******************************************************************************
  In addition to log likelihood, we also require a motif to actually be
  present a good portion of the time in an absolute sense across the genome.
  The coverage map is just a numerical map of whether or not to accept a
  putative motif as a real one (despite its log likelihood score.  A motif is
  considered "good" if it contains a 3-base subset of itself that is present
  in at least 20% of the total genes.  In the final stage of iterative start
  training, all motifs are labeled good.  0 = bad, 1 = good, 2 = good
  w/mismatch.
******************************************************************************/
void build_coverage_map(double real[4][4][4096], int good[4][4][4096], double
                        ng, int stage)
{
  int i, j, k, l, tmp, decomp[3];
  double thresh = 0.2;

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4096; k++)
      {
        good[i][j][k] = 0;
      }
    }
  }

  /* 3-base motifs */
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 64; j++)
    {
      if (real[0][i][j]/ng >= thresh)
      {
        for (k = 0; k < 4; k++)
        {
          good[0][k][j] = 1;
        }
      }
    }
  }

  /* 4-base motifs, must contain two valid 3-base motifs */
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 256; j++)
    {
      decomp[0] = (j&252)>>2;
      decomp[1] = j&63;
      if (good[0][i][decomp[0]] == 0 || good[0][i][decomp[1]] == 0)
      {
        continue;
      }
      good[1][i][j] = 1;
    }
  }

  /* 5-base motifs, interior mismatch allowed only if entire 5-base */
  /* motif represents 3 valid 3-base motifs (if mismatch converted) */
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 1024; j++)
    {
      decomp[0] = (j&1008)>>4;
      decomp[1] = (j&252)>>2;
      decomp[2] = j&63;
      if (good[0][i][decomp[0]] == 0 || good[0][i][decomp[1]] == 0 ||
          good[0][i][decomp[2]] == 0)
      {
        continue;
      }
      good[2][i][j] = 1;
      tmp = j;
      for (k = 0; k <= 16; k+= 16)
      {
        tmp = tmp ^ k;
        for (l = 0; l <= 32; l+= 32)
        {
          tmp = tmp ^ l;
          if (good[2][i][tmp] == 0)
          {
            good[2][i][tmp] = 2;
          }
        }
      }
    }
  }

  /* 6-base motifs, must contain two valid 5-base motifs */
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4096; j++)
    {
      decomp[0] = (j&4092)>>2;
      decomp[1] = j&1023;
      if (good[2][i][decomp[0]] == 0 || good[2][i][decomp[1]] == 0)
      {
        continue;
      }
      if (good[2][i][decomp[0]] == 1 && good[2][i][decomp[1]] == 1)
      {
        good[3][i][j] = 1;
      }
      else
      {
        good[3][i][j] = 2;
      }
    }
  }
}

/******************************************************************************
  Examines the results of the SD motif search to determine if this organism
  uses an SD motif or not.  Some motif of 3-6bp has to be good or we set
  uses_sd to 0, which will cause Prodigal to run the non-SD motif finder for
  starts.
******************************************************************************/
void determine_sd_usage(struct _training *tinf)
{
  tinf->uses_sd = 1;
  if (tinf->rbs_wt[0] >= 0.0)
  {
    tinf->uses_sd = 0;
  }
  if (tinf->rbs_wt[16] < 1.0 && tinf->rbs_wt[13] < 1.0 &&
      tinf->rbs_wt[15] < 1.0 && (tinf->rbs_wt[0] >= -0.5 ||
      (tinf->rbs_wt[22] < 2.0 && tinf->rbs_wt[24] < 2.0 &&
      tinf->rbs_wt[27] < 2.0)))
  {
    tinf->uses_sd = 0;
  }
}
