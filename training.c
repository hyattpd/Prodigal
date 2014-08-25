/******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2014 University of Tenum_nodesessee / UT-Battelle

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
int read_training_file(char *file_name, struct _training *train_data)
{
  size_t ret_val = 0;
  FILE *fh = NULL;

  fh = fopen(file_name, "rb");
  if (fh == NULL)
  {
    return -1;
  }
  ret_val = fread(train_data, sizeof(struct _training), 1, fh);
  fclose(fh);
  if (ret_val != 1)
  {
    return -1;
  }
  return 0;
}

/* Writes a training file to use for a later run of gene prediction */
int write_training_file(char *file_name, struct _training *train_data)
{
  size_t ret_val = 0;
  FILE *fh = NULL;

  fh = fopen(file_name, "wb");
  if (fh == NULL)
  {
    return -1;
  }
  ret_val = fwrite(train_data, sizeof(struct _training), 1, fh);
  if (ret_val != 1)
  {
    return -1;
  }
  return 0;
}

/* Build the training set, check its quality, and rebuild it with */
/* genetic code 4 if quality is low and in auto mode for genetic code */
void build_training_set_full(struct _node *nodes, struct _training *train_data,
                             struct _summary *statistics, unsigned char *seq,
                             unsigned char *rseq, unsigned char *useq,
                             int seq_length, int *num_nodes, int num_seq,
                             int genetic_code, int quiet)
{
  int train_qual = 0;
  char text[10000] = "";

  sprintf(text, "Building training set using genetic code %d...",
          train_data->trans_table);
  log_text(quiet, text);
  build_training_set(nodes, train_data, statistics, seq, rseq, useq,
                     seq_length, num_nodes, num_seq);
  log_text(quiet, "done!\n");

  /***********************************************************************
    Check average gene length to see if translation table looks good or
    if there is substantial gene decay.  If no genetic code was specified,
    try genetic code 4 to see if it solves the problem.
  ***********************************************************************/
  log_text(quiet, "Checking average training gene length...");
  train_qual = training_set_quality(statistics);
  /* Looks ok */
  if (train_qual == 0)
  {
    sprintf(text, "%.1f, looks ok.\n", statistics->avg_comp_gene_len);
    log_text(quiet, text);
  }
  /* Too many partial genes */
  else if (train_qual == 1)
  {
    log_text(quiet, "low but sequence is really drafty.\n");
    low_gene_len_warning(train_qual, statistics);
  }
  /* Poor average gene length */
  else if (train_qual == 2)
  {
    sprintf(text, "%.1f, too low.\n", statistics->avg_comp_gene_len);
    log_text(quiet, text);
    if (genetic_code == 0)
    {
      log_text(quiet, "Trying genetic code 4...");
      train_data->trans_table = 4;
      build_training_set(nodes, train_data, statistics, seq, rseq, useq,
                         seq_length, num_nodes, num_seq);
      train_qual = training_set_quality(statistics);
      if (train_qual < 2)
      {
        log_text(quiet, "looks good, using genetic code 4.\n");
      }
      else
      {
        log_text(quiet, "still bad, reverting to genetic code 11.\n");
        log_text(quiet, "Redoing genome with genetic code 11...");
        train_data->trans_table = 11;
        build_training_set(nodes, train_data, statistics, seq, rseq, useq,
                           seq_length, num_nodes, num_seq);
        log_text(quiet, "done.\n");
        low_gene_len_warning(train_qual, statistics);
      }
    }
    else
    {
      low_gene_len_warning(train_qual, statistics);
    }
  }
}

/* Build the training set using the supplied genetic code */
void build_training_set(struct _node *nodes, struct _training *train_data,
                        struct _summary *statistics, unsigned char *seq,
                        unsigned char *rseq, unsigned char *useq,
                        int seq_length, int *num_nodes, int num_seq)
{
  int initial_node = -1;
  int last_node = -1;

  /***********************************************************************
    Find all the potential starts and stops, sort them, and create a
    comprehensive list of nodes for dynamic programming.
  ***********************************************************************/
  if (*num_nodes > 0)
  {
    zero_nodes(nodes, *num_nodes);
  }
  *num_nodes = add_nodes(seq, rseq, useq, seq_length, nodes, 0, 0,
                         train_data->trans_table);
  qsort(nodes, *num_nodes, sizeof(struct _node), &compare_nodes);

  /***********************************************************************
    Scan all the ORFS looking for a potential GC bias in a particular
    codon position.  This information will be used to acquire a good
    initial set of genes.
  ***********************************************************************/
  record_gc_frame_bias(train_data, seq, seq_length, nodes, *num_nodes);

  /***********************************************************************
    Do an initial dynamic programming routine with just the GC frame
    bias used as a scoring function.  This will get an initial set of
    genes to train on.
  ***********************************************************************/
  record_overlapping_starts(nodes, *num_nodes, train_data->start_weight, 0);
  last_node = dynamic_programming(nodes, *num_nodes, train_data->bias,
                                  train_data->start_weight, 0);
  initial_node = find_first_node_from_last_node(nodes, last_node);

  /***********************************************************************
    Gather dicodon statistics for the training set.  Score the entire set
    of nodes.
  ***********************************************************************/
  calc_dicodon_gene(train_data, seq, rseq, seq_length, nodes, last_node);
  calc_coding_score(seq, rseq, seq_length, nodes, *num_nodes,
                    train_data->trans_table, train_data->gc,
                    train_data->gene_dc);

  /***********************************************************************
    Gather statistics about average gene length to see if the training
    set looks good.
  ***********************************************************************/
  calc_training_set_stats(nodes, initial_node, statistics, num_seq,
                          seq_length);
}

/* Records the GC frame bias from the node GC statistics */
void record_gc_frame_bias(struct _training *train_data, unsigned char *seq,
                          int seq_length, struct _node *nodes, int num_nodes)
{
  int i = 0;
  int len = 0;           /* Length of a gene */
  int *gc_frame = NULL;  /* Stores max frame for every position in sequence */
  double total = 0.0;    /* Sum of all the GC frame counts */

  /* Calculate highest GC frame for each position in sequence */
  gc_frame = calc_most_gc_frame(seq, seq_length);
  /* For each node, count all the occurrences of each of these frames */
  frame_score(gc_frame, nodes, num_nodes);
  free(gc_frame);

  /* Go through all the nodes and record all the bias info */
  /* Since every node is counted, sequence near the 3' edge gets */
  /* counted more often usually, and longer genes get weighted */
  /* more (similar to n-squared). */
  for (i = 0; i < 3; i++)
  {
    train_data->bias[i] = 0.0;
  }
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == START)
    {
      len = abs(nodes[i].stop_val-nodes[i].index)+1;
      train_data->bias[nodes[i].gc_bias] +=
        (nodes[i].gc_score[nodes[i].gc_bias]*len)/1000.0;
    }
  }
  total = train_data->bias[0] + train_data->bias[1] + train_data->bias[2];
  for (i = 0; i < 3; i++)
  {
    train_data->bias[i] *= (3.0/total);
  }
}

/******************************************************************************
  Simple routine that calculates the dicodon frequency in genes and in the
  background, and then stores the log likelihood of each 6-mer relative to the
  background.
******************************************************************************/
void calc_dicodon_gene(struct _training *train_data, unsigned char *seq,
                       unsigned char *rseq, int seq_length,
                       struct _node *nodes, int initial_node)
{
  int i = 0;
  int path = initial_node;         /* Dynamic programming traversal */
  int counts[4096] = {0};          /* Counts of each hexamer in real genes */
  int total = 0;                   /* Total of all hexamer counts */
  int left = -1;                   /* Left boundary of gene */
  int right = -1;                  /* Right boundary of gene */
  int in_gene = 0;                 /* -1 if not in gene, 1 if in gene */
  double prob[4096] = {0};         /* Counts divided by total */
  double background[4096] = {0};   /* Counts over the whole sequence */

  calc_mer_background(6, seq, rseq, seq_length, background);
  while (path != -1)
  {
    if (nodes[path].strand == -1 && nodes[path].type == START)
    {
      in_gene = -1;
      left = seq_length-nodes[path].index-1;
    }
    if (nodes[path].strand == 1 && nodes[path].type == STOP)
    {
      in_gene = 1;
      right = nodes[path].index+2;
    }
    if (in_gene == -1 && nodes[path].strand == -1 && nodes[path].type == STOP)
    {
      right = seq_length-nodes[path].index+1;
      for (i = left + SSTRUCT_SIZE; i < right-5; i+=3)
      {
        counts[mer_index(6, rseq, i)]++;
        total++;
      }
      in_gene = 0;
    }
    if (in_gene == 1 && nodes[path].strand == 1 && nodes[path].type == START)
    {
      left = nodes[path].index;
      for (i = left + SSTRUCT_SIZE; i < right-5; i+=3)
      {
        counts[mer_index(6, seq, i)]++;
        total++;
      }
      in_gene = 0;
    }
    path = nodes[path].trace_back;
  }
  for (i = 0; i < 4096; i++)
  {
    prob[i] = (counts[i]*1.0)/(total*1.0);
    if (prob[i] == 0 && background[i] != 0)
    {
      train_data->gene_dc[i] = -5.0;
    }
    else if (background[i] == 0)
    {
      train_data->gene_dc[i] = 0.0;
    }
    else
    {
      train_data->gene_dc[i] = log(prob[i]/background[i]);
    }
    if (train_data->gene_dc[i] > 5.0)
    {
      train_data->gene_dc[i] = 5.0;
    }
    if (train_data->gene_dc[i] < -5.0)
    {
      train_data->gene_dc[i] = -5.0;
    }
  }
}

/******************************************************************************
  Look at average complete gene length.  If it's above the minimum acceptable
  length, return 0.  If it's too small, see if this is due to the sequence
  being in tons of contigs.  If so, return 1.  If we still have no explanation
  for why it's low, return 2.
******************************************************************************/
int training_set_quality(struct _summary *genome_data)
{
  if (genome_data->avg_comp_gene_len > MIN_AVG_TRAIN_GENE_LEN)
  {
    return 0;
  }
  else if (genome_data->avg_contig_len < MIN_AVG_TRAIN_CTG_LEN ||
           genome_data->num_partial_genes > genome_data->num_complete_genes)
  {
    return 1;
  }
  else
  {
    return 2;
  }
}

/* Output a warning for low average gene length */
void low_gene_len_warning(int flag, struct _summary *genome_data)
{
  if (flag < 2)
  {
    fprintf(stderr, "\nWarning: Training sequence is highly fragmented.\n");
    fprintf(stderr, "You may get better results with the ");
    fprintf(stderr, "'-p anon' option.\n\n");
  }
  else
  {
    fprintf(stderr, "\nWarning: Average training gene length is");
    fprintf(stderr, " low (%.1f).\n", genome_data->avg_comp_gene_len);
    fprintf(stderr, "Double check translation table or check for");
    fprintf(stderr, " pseudogenes/gene decay.\n\n");
  }
}

/******************************************************************************
  Iterative Algorithm to train starts.  It begins with all the highest coding
  starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts moving
  starts around attempting to match these discoveries.  This start trainer is
  for Shine-Dalgarno motifs only.
******************************************************************************/
void train_starts_sd(unsigned char *seq, unsigned char *rseq, int seq_length,
                     struct _node *nodes, int num_nodes,
                     struct _training *train_data)
{
  int i = 0;               /* Loop variables */
  int j = 0;
  int frame = 0;           /* Reading frame */
  int best_rbs[3] = {0};   /* RBS of the best gene in each frame */
  int best_type[3] = {0};  /* Type of the best gene in each frame */
  int best_index[3] = {0}; /* Index of best gene in each frame */
  double best_score[3] = {0};  /* Best scoring node in each frame */
  int rbs_switch = 0;      /* Best motif between exact and mismatch */
  double denom = 0.0;      /* Denominator variable */
  double wt = 0.0;         /* Shorthand for start weight to save space */
  double rbg[28] = {0};    /* RBS background - i.e. all nodes */
  double rreal[28] = {0};  /* RBS foreground - i.e. only real genes */
  double tbg[3] = {0};     /* Type background - i.e. all nodes */
  double treal[3] = {0};   /* Type foreground - i.e. only real genes */
  double pcbg[48][4][4] = {{{0}}};  /* Pair composition background */
  double score_thresh = MIN_TRAIN_GENE_SCORE;  
                           /* Minimum score for genes to be "real" */

  wt = train_data->start_weight;
  zero_start_weights(train_data, 0);
  calc_type_background(nodes, num_nodes, tbg);
  calc_pair_comp_background(seq, rseq, seq_length, nodes, num_nodes, pcbg);

  /* Iterate SD_ITER times through the list of nodes                    */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few      */
  /* extra to be safe)                                                  */
  for (i = 0; i < SD_ITER; i++)
  {

    /* Recalculate the RBS motif background */
    for (j = 0; j < 28; j++)
    {
      rbg[j] = 0.0;
    }
    for (j = 0; j < num_nodes; j++)
    {
      if (nodes[j].type == STOP || nodes[j].edge == 1)
      {
        continue;
      }
      if (train_data->rbs_wt[nodes[j].rbs[0]] >
          train_data->rbs_wt[nodes[j].rbs[1]] + 1.0 ||
          nodes[j].rbs[1] == 0)
      {
        rbs_switch = nodes[j].rbs[0];
      }
      else if (train_data->rbs_wt[nodes[j].rbs[0]] <
               train_data->rbs_wt[nodes[j].rbs[1]] - 1.0 ||
               nodes[j].rbs[0] == 0)
      {
        rbs_switch = nodes[j].rbs[1];
      }
      else
      {
        rbs_switch = (int)dmax(nodes[j].rbs[0], nodes[j].rbs[1]);
      }
      rbg[rbs_switch] += 1.0;
    }
    denom = 0.0;
    for (j = 0; j < 28; j++)
    {
      denom += rbg[j];
    }
    for (j = 0; j < 28; j++)
    {
      rbg[j] /= denom;
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
      best_score[j] = 0.0;
      best_index[j] = -1;
      best_rbs[j] = 0;
      best_type[j] = 0;
    }
    for (j = 0; j < num_nodes; j++)
    {
      if (nodes[j].type == START && nodes[j].edge == 1)
      {
        continue;
      }
      frame = (nodes[j].index)%3;
      if (nodes[j].type == STOP && nodes[j].strand == 1)
      {
        if (best_score[frame] >= score_thresh &&
            nodes[best_index[frame]].index%3 == frame)
        {
          rreal[best_rbs[frame]] += 1.0;
          treal[best_type[frame]] += 1.0;
          if (i == SD_ITER-1)
          {
            count_upstream_composition(seq, seq_length, 1,
                                       nodes[best_index[frame]].index,
                                       train_data);
            count_pair_composition(seq, seq_length, 1,
                                   nodes[best_index[frame]].index,
                                   train_data->pair_comp);
          }
        }
        best_score[frame] = 0.0;
        best_index[frame] = -1;
        best_rbs[frame] = 0;
        best_type[frame] = 0;
      }
      else if (nodes[j].strand == 1)
      {
        if (train_data->rbs_wt[nodes[j].rbs[0]] >
            train_data->rbs_wt[nodes[j].rbs[1]] + 1.0 ||
            nodes[j].rbs[1] == 0)
        {
          rbs_switch = nodes[j].rbs[0];
        }
        else if (train_data->rbs_wt[nodes[j].rbs[0]] <
                 train_data->rbs_wt[nodes[j].rbs[1]] - 1.0 ||
                 nodes[j].rbs[0] == 0)
        {
          rbs_switch = nodes[j].rbs[1];
        }
        else
        {
          rbs_switch = (int)dmax(nodes[j].rbs[0], nodes[j].rbs[1]);
        }
        if (nodes[j].cscore + wt*train_data->rbs_wt[rbs_switch] +
            wt*train_data->type_wt[nodes[j].subtype] >= best_score[frame])
        {
          best_score[frame] = nodes[j].cscore + wt*train_data->rbs_wt[rbs_switch];
          best_score[frame] += wt*train_data->type_wt[nodes[j].subtype];
          best_index[frame] = j;
          best_type[frame] = nodes[j].subtype;
          best_rbs[frame] = rbs_switch;
        }
      }
    }

    /* Reverse strand pass */
    for (j = 0; j < 3; j++)
    {
      best_score[j] = 0.0;
      best_index[j] = -1;
      best_rbs[j] = 0;
      best_type[j] = 0;
    }
    for (j = num_nodes-1; j >= 0; j--)
    {
      if (nodes[j].type == START && nodes[j].edge == 1)
      {
        continue;
      }
      frame = (nodes[j].index)%3;
      if (nodes[j].type == STOP && nodes[j].strand == -1)
      {
        if (best_score[frame] >= score_thresh && nodes[best_index[frame]].index%3 == frame)
        {
          rreal[best_rbs[frame]] += 1.0;
          treal[best_type[frame]] += 1.0;
          if (i == SD_ITER-1)
          {
            count_upstream_composition(rseq, seq_length, -1,
                                       nodes[best_index[frame]].index, train_data);
            count_pair_composition(rseq, seq_length, -1,
                                   nodes[best_index[frame]].index,
                                   train_data->pair_comp);
          }
        }
        best_score[frame] = 0.0;
        best_index[frame] = -1;
        best_rbs[frame] = 0;
        best_type[frame] = 0;
      }
      else if (nodes[j].strand == -1)
      {
        if (train_data->rbs_wt[nodes[j].rbs[0]] >
            train_data->rbs_wt[nodes[j].rbs[1]] + 1.0 ||
            nodes[j].rbs[1] == 0)
        {
          rbs_switch = nodes[j].rbs[0];
        }
        else if (train_data->rbs_wt[nodes[j].rbs[0]] <
                 train_data->rbs_wt[nodes[j].rbs[1]] - 1.0 ||
                 nodes[j].rbs[0] == 0)
        {
          rbs_switch = nodes[j].rbs[1];
        }
        else
        {
          rbs_switch = (int)dmax(nodes[j].rbs[0], nodes[j].rbs[1]);
        }
        if (nodes[j].cscore + wt*train_data->rbs_wt[rbs_switch] +
            wt*train_data->type_wt[nodes[j].subtype] >= best_score[frame])
        {
          best_score[frame] = nodes[j].cscore + wt*train_data->rbs_wt[rbs_switch];
          best_score[frame] += wt*train_data->type_wt[nodes[j].subtype];
          best_index[frame] = j;
          best_type[frame] = nodes[j].subtype;
          best_rbs[frame] = rbs_switch;
        }
      }
    }

    denom = 0.0;
    for (j = 0; j < 28; j++)
    {
      denom += rreal[j];
    }
    if (denom == 0.0)
    {
      for (j = 0; j < 28; j++)
      {
        train_data->rbs_wt[j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 28; j++)
      {
        rreal[j] /= denom;
        if (rbg[j] != 0)
        {
          train_data->rbs_wt[j] = log(rreal[j]/rbg[j]);
        }
        else
        {
          train_data->rbs_wt[j] = -4.0;
        }
        if (train_data->rbs_wt[j] > 4.0)
        {
          train_data->rbs_wt[j] = 4.0;
        }
        if (train_data->rbs_wt[j] < -4.0)
        {
          train_data->rbs_wt[j] = -4.0;
        }
      }
    }
    denom = 0.0;
    for (j = 0; j < 3; j++)
    {
      denom += treal[j];
    }
    if (denom == 0.0)
    {
      for (j = 0; j < 3; j++)
      {
        train_data->type_wt[j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 3; j++)
      {
        treal[j] /= denom;
        if (tbg[j] != 0)
        {
          train_data->type_wt[j] = log(treal[j]/tbg[j]);
        }
        else
        {
          train_data->type_wt[j] = -4.0;
        }
        if (train_data->type_wt[j] > 4.0)
        {
          train_data->type_wt[j] = 4.0;
        }
        if (train_data->type_wt[j] < -4.0)
        {
          train_data->type_wt[j] = -4.0;
        }
      }
    }
    if (denom <= (double)num_nodes/SMALL_TRAIN_SET_NODES)
    {
      score_thresh /= 2.0;
    }
  }

  /* Convert upstream base composition to a log score */
  for (i = 0; i < 32; i++)
  {
    denom = 0.0;
    for (j = 0; j < 4; j++)
    {
      denom += train_data->ups_comp[i][j];
    }
    if (denom == 0.0)
    {
      for (j = 0; j < 4; j++)
      {
        train_data->ups_comp[i][j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 4; j++)
      {
        train_data->ups_comp[i][j] /= denom;
        if (train_data->gc > 0.1 && train_data->gc < 0.9)
        {
          if (j == 0 || j == 3)
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/(1.0-train_data->gc));
          }
          else
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/train_data->gc);
          }
        }
        else if (train_data->gc <= 0.1)
        {
          if (j == 0 || j == 3)
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.90);
          }
          else
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.10);
          }
        }
        else
        {
          if (j == 0 || j == 3)
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.10);
          }
          else
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.90);
          }
        }
        if (train_data->ups_comp[i][j] > 4.0)
        {
          train_data->ups_comp[i][j] = 4.0;
        }
        if (train_data->ups_comp[i][j] < -4.0)
        {
          train_data->ups_comp[i][j] = -4.0;
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
void train_starts_nonsd(unsigned char *seq, unsigned char *rseq,
                        int seq_length, struct _node *nodes, int num_nodes,
                        struct _training *train_data)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int fr = 0;              /* Reading frame */
  int stage = 0;           /* Motifs go through 3 stages of counting */
  int best_index[3] = {0}; /* Index of best gene in each frame */
  int good_motif[4][4][4096] = {{{0}}};  /* RBS post-coverage-check */
  double denom = 0.0;      /* Denominator variable */
  double wt = 0.0;         /* Shorthand for start weight to save space */
  double mbg[4][4][4096] = {{{0}}};    /* RBS background - i.e. all nodes */
  double mreal[4][4][4096] = {{{0}}};  /* RBS foreground - i.e. real genes */
  double best_score[3] = {0};    /* Best scoring node in each frame */
  double tbg[3] = {0};     /* Type background - i.e. all nodes */
  double treal[3] = {0};   /* Type foreground - i.e. only real genes */
  double pcbg[48][4][4] = {{{0}}};  /* Pair composition background */
  double zbg = 0.0;        /* Background for zero motif */
  double zreal = 0.0;      /* Zero motif weight for real genes */
  double ngenes = 0.0;     /* Gene counter, double in case of weighting */
  double score_thresh = MIN_TRAIN_GENE_SCORE;
                           /* Minimum score for genes to be "real" */

  wt = train_data->start_weight;
  zero_start_weights(train_data, 1);
  calc_type_background(nodes, num_nodes, tbg);
  calc_pair_comp_background(seq, rseq, seq_length, nodes, num_nodes, pcbg);

  /* Iterate NONSD_ITER times through the list of nodes                 */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few      */
  /* extra to be safe)                                                  */
  for (i = 0; i < NONSD_ITER; i++)
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
    for (j = 0; j < num_nodes; j++)
    {
      if (nodes[j].type == STOP || nodes[j].edge == 1)
      {
        continue;
      }
      find_best_nonsd_motif(train_data, seq, rseq, seq_length, &nodes[j],
                            stage);
      update_nonsd_motif_counts(mbg, &zbg, seq, rseq, seq_length, &(nodes[j]),
                                stage);
    }
    denom = 0.0;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          denom += mbg[j][k][l];
        }
      }
    }
    denom += zbg;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          mbg[j][k][l] /= denom;
        }
      }
    }
    zbg /= denom;

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
      best_score[j] = 0.0;
      best_index[j] = -1;
    }
    for (j = 0; j < num_nodes; j++)
    {
      if (nodes[j].type == START && nodes[j].edge == 1)
      {
        continue;
      }
      fr = (nodes[j].index)%3;
      if (nodes[j].type == STOP && nodes[j].strand == 1)
      {
        if (best_score[fr] >= score_thresh)
        {
          ngenes += 1.0;
          treal[nodes[best_index[fr]].subtype] += 1.0;
          update_nonsd_motif_counts(mreal, &zreal, seq, rseq, seq_length,
                                    &(nodes[best_index[fr]]), stage);
          if (i == NONSD_ITER-1)
          {
            count_upstream_composition(seq, seq_length, 1,
                                       nodes[best_index[fr]].index, train_data);
            count_pair_composition(seq, seq_length, 1,
                                   nodes[best_index[fr]].index,
                                   train_data->pair_comp);
          }
        }
        best_score[fr] = 0.0;
        best_index[fr] = -1;
      }
      else if (nodes[j].strand == 1)
      {
        if (nodes[j].cscore + wt*nodes[j].mot.score +
            wt*train_data->type_wt[nodes[j].subtype]
            >= best_score[fr])
        {
          best_score[fr] = nodes[j].cscore + wt*nodes[j].mot.score;
          best_score[fr] += wt*train_data->type_wt[nodes[j].subtype];
          best_index[fr] = j;
        }
      }
    }

    /* Reverse strand pass */
    for (j = 0; j < 3; j++)
    {
      best_score[j] = 0.0;
      best_index[j] = -1;
    }
    for (j = num_nodes-1; j >= 0; j--)
    {
      if (nodes[j].type == START && nodes[j].edge == 1)
      {
        continue;
      }
      fr = (nodes[j].index)%3;
      if (nodes[j].type == STOP && nodes[j].strand == -1)
      {
        if (best_score[fr] >= score_thresh)
        {
          ngenes += 1.0;
          treal[nodes[best_index[fr]].subtype] += 1.0;
          update_nonsd_motif_counts(mreal, &zreal, seq, rseq, seq_length,
                                    &(nodes[best_index[fr]]), stage);
          if (i == NONSD_ITER-1)
          {
            count_upstream_composition(rseq, seq_length, -1,
                                       nodes[best_index[fr]].index, train_data);
            count_pair_composition(rseq, seq_length, -1,
                                   nodes[best_index[fr]].index,
                                   train_data->pair_comp);
          }
        }
        best_score[fr] = 0.0;
        best_index[fr] = -1;
      }
      else if (nodes[j].strand == -1)
      {
        if (nodes[j].cscore + wt*nodes[j].mot.score +
            wt*train_data->type_wt[nodes[j].subtype]
            >= best_score[fr])
        {
          best_score[fr] = nodes[j].cscore + wt*nodes[j].mot.score;
          best_score[fr] += wt*train_data->type_wt[nodes[j].subtype];
          best_index[fr] = j;
        }
      }
    }

    /* Update the log likelihood weights for type and RBS motifs */
    if (stage < 2)
    {
      label_good_nonsd_motifs(mreal, good_motif, ngenes);
    }
    denom = 0.0;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        for (l = 0; l < 4096; l++)
        {
          denom += mreal[j][k][l];
        }
      }
    }
    denom += zreal;
    if (denom == 0.0)
    {
      for (j = 0; j < 4; j++)
      {
        for (k = 0; k < 4; k++)
        {
          for (l = 0; l < 4096; l++)
          {
            train_data->mot_wt[j][k][l] = 0.0;
          }
        }
      }
      train_data->no_mot = 0.0;
    }
    else
    {
      for (j = 0; j < 4; j++)
      {
        for (k = 0; k < 4; k++)
        {
          for (l = 0; l < 4096; l++)
          {
            if (good_motif[j][k][l] == 0)
            {
              zreal += mreal[j][k][l];
              zbg += mreal[j][k][l];
              mreal[j][k][l] = 0.0;
              mbg[j][k][l] = 0.0;
            }
            mreal[j][k][l] /= denom;
            if (mbg[j][k][l] != 0)
            {
              train_data->mot_wt[j][k][l] = log(mreal[j][k][l]/mbg[j][k][l]);
            }
            else
            {
              train_data->mot_wt[j][k][l] = -4.0;
            }
            if (train_data->mot_wt[j][k][l] > 4.0)
            {
              train_data->mot_wt[j][k][l] = 4.0;
            }
            if (train_data->mot_wt[j][k][l] < -4.0)
            {
              train_data->mot_wt[j][k][l] = -4.0;
            }
          }
        }
      }
    }
    zreal /= denom;
    if (zbg != 0)
    {
      train_data->no_mot = log(zreal/zbg);
    }
    else
    {
      train_data->no_mot = -4.0;
    }
    if (train_data->no_mot > 4.0)
    {
      train_data->no_mot = 4.0;
    }
    if (train_data->no_mot < -4.0)
    {
      train_data->no_mot = -4.0;
    }
    denom = 0.0;
    for (j = 0; j < 3; j++)
    {
      denom += treal[j];
    }
    if (denom == 0.0)
    {
      for (j = 0; j < 3; j++)
      {
        train_data->type_wt[j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 3; j++)
      {
        treal[j] /= denom;
        if (tbg[j] != 0)
        {
          train_data->type_wt[j] = log(treal[j]/tbg[j]);
        }
        else
        {
          train_data->type_wt[j] = -4.0;
        }
        if (train_data->type_wt[j] > 4.0)
        {
          train_data->type_wt[j] = 4.0;
        }
        if (train_data->type_wt[j] < -4.0)
        {
          train_data->type_wt[j] = -4.0;
        }
      }
    }
    if (denom <= (double)num_nodes/SMALL_TRAIN_SET_NODES)
    {
      score_thresh /= 2.0;
    }
  }
  /* Convert upstream base composition to a log score */
  for (i = 0; i < 32; i++)
  {
    denom = 0.0;
    for (j = 0; j < 4; j++)
    {
      denom += train_data->ups_comp[i][j];
    }
    if (denom == 0.0)
    {
      for (j = 0; j < 4; j++)
      {
        train_data->ups_comp[i][j] = 0.0;
      }
    }
    else
    {
      for (j = 0; j < 4; j++)
      {
        train_data->ups_comp[i][j] /= denom;
        if (train_data->gc > 0.1 && train_data->gc < 0.9)
        {
          if (j == 0 || j == 3)
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/(1.0-train_data->gc));
          }
          else
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/train_data->gc);
          }
        }
        else if (train_data->gc <= 0.1)
        {
          if (j == 0 || j == 3)
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.90);
          }
          else
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.10);
          }
        }
        else
        {
          if (j == 0 || j == 3)
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.10);
          }
          else
          {
            train_data->ups_comp[i][j] =
              log(train_data->ups_comp[i][j]*2.0/0.90);
          }
        }
        if (train_data->ups_comp[i][j] > 4.0)
        {
          train_data->ups_comp[i][j] = 4.0;
        }
        if (train_data->ups_comp[i][j] < -4.0)
        {
          train_data->ups_comp[i][j] = -4.0;
        }
      }
    }
  }
}

/******************************************************************************
  For a given start, record the base composition of the 24 bases downstream
  and upstream of the start site, as well as the bases they could potentially
  pair with to produce secondary structure.  Returns a 1 on success, 0 on
  failure (return value is meant to be used as a count if succeeded).
******************************************************************************/
int count_pair_composition(unsigned char *seq, int seq_length, int strand,
                           int pos, double pcdata[48][4][4])
{
  int i = 0;
  int j = 0;
  int start = 0;
  int count[4] = {0};         /* Simple count of A, C, T, G */
  int index = 0;              /* Index to hold mer index of base */
  int prohibit = 0;           /* Index of pairs that can't bond */

  if (strand == 1)
  {
    start = pos;
  }
  else
  {
    start = seq_length-1-pos;
  }

  if (start - SSTRUCT_SIZE < 0 || start + SSTRUCT_SIZE > seq_length)
  {
    return 0;
  }

  /* Gather simple count of ACTG in total window */
  for (i = -SSTRUCT_SIZE; i < SSTRUCT_SIZE; i++)
  {
    count[mer_index(1, seq, start+i)]++;
  }

  /* Now do the pair count */
  for (i = -SSTRUCT_SIZE; i < SSTRUCT_SIZE; i++)
  {
    index = mer_index(1, seq, start+i);
    for (j = 0; j < 4; j++)
    {
      pcdata[i+SSTRUCT_SIZE][index][j] += count[j];
    }
    for (j = -3; j < 3; j++)
    {
      if (i + j < -SSTRUCT_SIZE || i + j >= SSTRUCT_SIZE)
      {
        continue;
      }
      prohibit = mer_index(1, seq, start+i+j);
      pcdata[i+SSTRUCT_SIZE][index][prohibit]--;
    }
  }

  return 1;
}

/******************************************************************************
  For a given start, record the base composition of the upstream region at
  positions -1 and -2 and -15 to -44.  This will be used to supplement the
  SD (or other) motif finder with additional information.
******************************************************************************/
void count_upstream_composition(unsigned char *seq, int seq_length, int strand,
                                int pos, struct _training *train_data)
{
  int i = 0;
  int start = 0;
  int count = 0;
  if (strand == 1)
  {
    start = pos;
  }
  else
  {
    start = seq_length-1-pos;
  }

  if (start - 45 < 0)
  {
    return;
  }

  /* At this stage, ups_comp just stores simple counts at each position */
  for (i = 1; i < 45; i++)
  {
    if (i > 2 && i < 15)
    {
      continue;
    }
    train_data->ups_comp[count][mer_index(1, seq, start-i)]++;
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
void update_nonsd_motif_counts(double motifs[4][4][4096], double *zero,
                               unsigned char *seq, unsigned char *rseq,
                               int seq_length, struct _node *nodes, int stage)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int start = 0;               /* Position of node in working sequence */
  int spacer_index = 0;        /* Spacer index in weight array */
  unsigned char *wseq = NULL;  /* Working sequence - seq or rseq */
  struct _motif *mot = &(nodes->mot);   /* Working motif structure */

  if (nodes->type == STOP || nodes->edge == 1)
  {
    return;
  }
  if (mot->len == 0)
  {
    *zero += 1.0;
    return;
  }

  if (nodes->strand == 1)
  {
    wseq = seq;
    start = nodes->index;
  }
  else
  {
    wseq = rseq;
    start = seq_length-1-nodes->index;
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
          spacer_index = 3;
        }
        else if (j <= start-14-i)
        {
          spacer_index = 2;
        }
        else if (j >= start-7-i)
        {
          spacer_index = 1;
        }
        else
        {
          spacer_index = 0;
        }
        for (k = 0; k < 4; k++)
        {
          motifs[i][k][mer_index(i+3, wseq, j)] += 1.0;
        }
      }
    }
  }
  /* Stage 1:  Count only the best motif, but also count  */
  /* all its sub-motifs.                                  */
  else if (stage == 1)
  {
    motifs[mot->len-3][mot->spacer_index][mot->index] += 1.0;
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
          spacer_index = 3;
        }
        else if (j <= start-14-i)
        {
          spacer_index = 2;
        }
        else if (j >= start-7-i)
        {
          spacer_index = 1;
        }
        else
        {
          spacer_index = 0;
        }
        motifs[i][spacer_index][mer_index(i+3, wseq, j)] += 1.0;
      }
    }
  }
  /* Stage 2:  Only count the highest scoring motif. */
  else if (stage == 2)
  {
    motifs[mot->len-3][mot->spacer_index][mot->index] += 1.0;
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
void label_good_nonsd_motifs(double real[4][4][4096], int good[4][4][4096],
                             double num_genes)
{
  int i = 0;                /* Loop variables */
  int j = 0;
  int k = 0;
  int l = 0;
  int tmp = 0;              /* Temporary variable to hold indices */
  int decomp[3] = {0};      /* Used to extract bits from the motif */
  double cvg_thresh = 0.2;  /* Coverage threshold */

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
      if (real[0][i][j]/num_genes >= cvg_thresh)
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
void determine_sd_usage(struct _training *train_data)
{
  train_data->uses_sd = 1;
  if (train_data->rbs_wt[0] >= 0.0)
  {
    train_data->uses_sd = 0;
  }
  if (train_data->rbs_wt[16] < 1.0 && train_data->rbs_wt[13] < 1.0 &&
      train_data->rbs_wt[15] < 1.0 && (train_data->rbs_wt[0] >= -0.5 ||
      (train_data->rbs_wt[22] < 2.0 && train_data->rbs_wt[24] < 2.0 &&
      train_data->rbs_wt[27] < 2.0)))
  {
    train_data->uses_sd = 0;
  }
}

/* Zero out start weights, 'which' is SD vs. nonSD */
void zero_start_weights(struct _training *train_data, int which)
{
  int i = 0;
  int j = 0;
  int k = 0;

  for (i = 0; j < 3; j++)
  {
    train_data->type_wt[j] = 0.0;
  }
  for (i = 0; i < 48; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        train_data->pair_comp[i][j][k] = 0.0;
      }
    }
  }

  if (which == 0)
  {
    for (i = 0; i < 28; i++)
    {
      train_data->rbs_wt[j] = 0.0;
    }
  }
  else if (which == 1)
  {
    train_data->no_mot = 0.0;
    for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
      {
        for (k = 0; k < 4096; k++)
        {
          train_data->mot_wt[i][j][k] = 0.0;
        }
      }   
    }
  }
}

/* Calculate ATG/GTG/TTG for every single node and build background */
void calc_type_background(struct _node *nodes, int num_nodes, double *tbg)
{
  int i = 0;
  int denom = 0;

  /* Build the background of random types */
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP)
    {
      continue;
    }
    tbg[nodes[i].subtype] += 1.0;
  }
  denom = 0.0;
  for (i = 0; i < 3; i++)
  {
    denom += tbg[i];
  }
  for (i = 0; i < 3; i++)
  {
    tbg[i] /= denom;
  }
}

/* Calculate the pair composition background */
void calc_pair_comp_background(unsigned char *seq, unsigned char *rseq,
                               int seq_length, struct _node *nodes,
                               int num_nodes, double pcbg[48][4][4])
{
  int i = 0;
  int j = 0;
  int k = 0;
  double denom = 0;

  /* Build the background of pair composition */
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP)
    {
      continue;
    }
    if (nodes[i].strand == 1)
    {
      denom += (double)count_pair_composition(seq, seq_length, 1,
                                              nodes[i].index, pcbg);
    }
    else if (nodes[i].strand == -1)
    {
      denom += (double)count_pair_composition(rseq, seq_length, -1,
                                              nodes[i].index, pcbg);
    }
  }

  for (i = 0; i < 48; i++)
  {
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        pcbg[i][j][k] /= denom;
      }
    }
  }

  for (i = 0; i < 48; i++)
  {
    denom = 0.0;
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        denom += pcbg[i][j][k];
      }
    }
    for (j = 0; j < 4; j++)
    {
      for (k = 0; k < 4; k++)
      {
        pcbg[i][j][k] /= denom;
/*        printf("%d\t%d\t%d\t%.6f\n", i, j, k, pcbg[i][j][k]); */
      }
    }
  }
}
