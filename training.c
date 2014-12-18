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
  zero_nodes(nodes, *num_nodes);
  *num_nodes = add_nodes(seq, rseq, useq, seq_length, nodes, 0, 0,
                         train_data->trans_table);
  qsort(nodes, *num_nodes, sizeof(struct _node), &compare_nodes);

  /* Get base probability of stop and start codons */
  calc_start_and_stop_probs(seq, rseq, seq_length, train_data);

  /***********************************************************************
    Scan all the ORFS looking for a potential GC bias in a particular
    codon position.  This information will be used to acquire a good
    initial set of genes.
  ***********************************************************************/
  frame_plot_score(seq, seq_length, nodes, *num_nodes, train_data);

  /***********************************************************************
    Do an initial dynamic programming routine with just the GC frame
    bias used as a scoring function.  This will get an initial set of
    genes to train on.
  ***********************************************************************/
  record_overlapping_starts(nodes, *num_nodes, train_data->start_weight, 0);
  last_node = dynamic_programming(nodes, *num_nodes, train_data->start_weight,
                                  0);
  initial_node = find_first_node_from_last_node(nodes, last_node);

  /***********************************************************************
    Gather dicodon statistics for the training set.  Score the entire set
    of nodes.
  ***********************************************************************/
  calc_dicodon_gene(train_data, seq, rseq, seq_length, nodes, last_node);
  calc_coding_score(seq, rseq, seq_length, nodes, *num_nodes, train_data);

  /***********************************************************************
    Gather statistics about average gene length to see if the training
    set looks good.
  ***********************************************************************/
  calc_training_set_stats(nodes, initial_node, statistics, num_seq,
                          seq_length);
}

/* Calculates base probabilities a codon is a start/stop */
void calc_start_and_stop_probs(unsigned char *seq, unsigned char *rseq,
                               int seq_length, struct _training *train_data)
{
  int i = 0;
  int trans_table = train_data->trans_table;
  double total_ctr = 0.0;
  double start_ctr = 0.0;
  double stop_ctr = 0.0;

  for (i = 0; i <= seq_length - 3; i++)
  {
    total_ctr += 2.0;
    if (get_start_type(seq, i, trans_table) != -1)
    {
      start_ctr += 1.0;
    }
    if (get_start_type(rseq, i, trans_table) != -1)
    {
      start_ctr += 1.0;
    }
    if (get_stop_type(seq, i, trans_table) != -1)
    {
      stop_ctr += 1.0;
    }
    if (get_stop_type(rseq, i, trans_table) != -1)
    {
      stop_ctr += 1.0;
    }
  }
  train_data->prob_start = start_ctr / total_ctr;
  train_data->prob_stop = stop_ctr / total_ctr;
}

/* Records the GC frame bias from the node GC statistics */
void frame_plot_score(unsigned char *seq, int seq_length, struct _node *nodes,
                      int num_nodes, struct _training *train_data)
{
  int i = 0;
  int j = 0;
  int *gc_frame = NULL;  /* Stores max frame for every position in sequence */
  double bias[3] = {0};  /* Frame bias log weights */
  int counter[3][3] = {{0}};     /* Counts for winning frame in each frame */
  int last[3] = {0};             /* Last stop codon we saw in each frame */
  int frame_mod = 0;             /* Relative-to-absolute frame modifier */
  int frame = 0;                 /* Frame of the current node */
  int best_frame = 0;            /* The final winning frame for a node */
  double gene_len = 0.0;         /* Length of the gene */
  double base_prob = 0.0;        /* Base probability a node is a gene */

  base_prob = 1.0 / (1100 * (train_data->prob_stop) * 6.0);

  /* Calculate highest GC frame for each position in sequence */
  gc_frame = calc_most_gc_frame(seq, seq_length);

  /* For each node, count all the occurrences of each of these frames */

  /* Go through all the nodes and record all the bias info */
  /* Since every node is counted, sequence near the 3' edge gets */
  /* counted more often usually, and longer genes get weighted */
  /* more (similar to n-squared). */
  /* Forward strand nodes */
  for (i = num_nodes-1; i >= 0; i--)
  {
    frame = (nodes[i].index)%3;
    frame_mod = 3 - frame;
    if (nodes[i].strand == 1 && nodes[i].type == STOP)
    {
      for (j = 0; j < 3; j++)
      {
        counter[frame][j] = 0;
      }
      last[frame] = nodes[i].index;
      counter[frame][(gc_frame[nodes[i].index] + frame_mod)%3] = 1;
    }
    else if (nodes[i].strand == 1)
    {
      for (j = last[frame]-3; j >= nodes[i].index; j-=3)
      {
        counter[frame][(gc_frame[j] + frame_mod)%3]++;
      }
      for (j = 0; j < 3; j++)
      {
        nodes[i].gc_frame[j] += counter[frame][j];
      }
      best_frame = max_frame(counter[frame][0], counter[frame][1],
                             counter[frame][2]);
      bias[best_frame] += counter[frame][best_frame]/1000.0;
      last[frame] = nodes[i].index;
    }
  }
  /* Reverse Strand Nodes */
  for (i = 0; i < num_nodes; i++)
  {
    frame = (nodes[i].index)%3;
    frame_mod = frame;
    if (nodes[i].strand == -1 && nodes[i].type == STOP)
    {
      for (j = 0; j < 3; j++)
      {
        counter[frame][j] = 0;
      }
      last[frame] = nodes[i].index;
      counter[frame][((3-gc_frame[nodes[i].index]) + frame_mod)%3] = 1;
    }
    else if (nodes[i].strand == -1)
    {
      for (j = last[frame]+3; j <= nodes[i].index; j+=3)
      {
        counter[frame][((3-gc_frame[j]) + frame_mod)%3]++;
      }
      for (j = 0; j < 3; j++)
      {
        nodes[i].gc_frame[j] += counter[frame][j];
      }
      best_frame = max_frame(counter[frame][0], counter[frame][1],
                             counter[frame][2]);
      bias[best_frame] += counter[frame][best_frame]/1000.0;
      last[frame] = nodes[i].index;
    }
  }
  free(gc_frame);

  /* Further normalization to reduce the log bias terms to sum to 1 */
  normalize_array(bias, 3);

  /* Now create a coding score for each node.  Score uses length */
  /* and the frame bias.  Effective gene length is modified by the */
  /* frame bias to be longer or shorter. */
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP)
    {
      continue;
    }
    gene_len = 3.0 * (bias[0]*nodes[i].gc_frame[0] +
                      bias[1]*nodes[i].gc_frame[1] +
                      bias[2]*nodes[i].gc_frame[2]);
    nodes[i].cscore = (1-gene_len) * log(1.0-train_data->prob_stop);
    nodes[i].cscore -= log(train_data->prob_start);
    nodes[i].cscore += log(base_prob);
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
  int j = 0;
  int path = initial_node;         /* Dynamic programming traversal */
  int left = -1;                   /* Left boundary of gene */
  int right = -1;                  /* Right boundary of gene */
  int in_gene = 0;                 /* -1 if not in gene, 1 if in gene */
  int index = 0;                   /* Index in 4096-length array */
  double prob = 0.0;               /* Probability a gene is real */
  double sum_real = 0.0;           /* Total of all real hexamer counts */
  double sum_bg = 0.0;             /* Total of all background hexamer counts */
  double counts[4096] = {0};       /* Counts within a single gene */
  double background[4096] = {0};   /* Counts over the whole sequence */
  double stop_bg[4] = {0};         /* Background for stop codons */
  double stop_real[4] = {0};       /* Stop codon counts in real genes */
 char qt[10];
 
  /* Count words across whole sequence */
  get_word_counts(6, seq, rseq, seq_length, background);
  calc_stop_background(seq, rseq, seq_length, train_data->trans_table,
                       stop_bg);

  /* Count words in training set genes */
  while (path != -1)
  {
    if (nodes[path].strand == -1 && nodes[path].type == START)
    {
      in_gene = -1;
      left = seq_length-nodes[path].index-1;
      prob = calculate_confidence(nodes[path].cscore, 1.0)/100.0;
    }
    if (nodes[path].strand == 1 && nodes[path].type == STOP)
    {
      in_gene = 1;
      right = nodes[path].index+2;
      stop_real[nodes[path].subtype] += 1.0;
    }
    if (in_gene == -1 && nodes[path].strand == -1 && nodes[path].type == STOP)
    {
      stop_real[nodes[path].subtype] += 1.0;
      right = seq_length-nodes[path].index+1;
      for (i = left; i < right - 5; i+=3)
      {
        index = mer_index(6, rseq, i);
        counts[index] += prob;
        background[index] -= prob;
      }
      in_gene = 0;
    }
    if (in_gene == 1 && nodes[path].strand == 1 && nodes[path].type == START)
    {
      prob = calculate_confidence(nodes[path].cscore, 1.0)/100.0;
      left = nodes[path].index;
      for (i = left; i < right - 8; i+=3)
      {
        index = mer_index(6, seq, i);
        counts[index] += prob;
        background[index] -= prob;
      }
      in_gene = 0;
    }
    path = nodes[path].trace_back;
  }
  /* Convert counts to normalized frequency */
/*
  for (i = 0; i < 4096; i++)
  {
    sum_real += counts[i];
    sum_bg += background[i];
  }
  for (i = 0; i < 4096; i++)
  {
    counts[i] /= sum_real;
    background[i] /= sum_bg;
    train_data->gene_dc[i] = log(counts[i]/background[i]);
  }
*/
printf("stop weights: %.4f\t%.4f\t%.4f\t%.4f\n", stop_bg[0], stop_bg[1], stop_bg[2], stop_bg[3]);
  normalize_array(stop_real, 4);
printf("stop weights: %.4f\t%.4f\t%.4f\t%.4f\n", stop_real[0], stop_real[1], stop_real[2], stop_real[3]);
  create_log_score(train_data->stop_wt, stop_real, stop_bg, 4);
printf("stop weights: %.4f\t%.4f\t%.4f\t%.4f\n", train_data->stop_wt[0], train_data->stop_wt[1], train_data->stop_wt[2], train_data->stop_wt[3]);
  for (i = 0; i < 64; i++)
  {
    sum_real = 0.0;
    sum_bg = 0.0;
    for (j = 0; j < 64; j++)
    {
      index = j*64 + i;
      if (counts[index] > 0 && background[index] > 0)
      {
        sum_real += counts[index];
        sum_bg += background[index];
      }
    }
    for (j = 0; j < 64; j++)
    {
      index = j*64 + i;
      if (counts[index] > 0 && background[index] > 0)
      {
        counts[index] /= sum_real;
        background[index] /= sum_bg;
        train_data->gene_dc[index] = log(counts[index]/background[index]);
      }
      else
      {
        train_data->gene_dc[index] = -4.0;
      }
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
  int k = 0;
  int num_genes = 0;       /* Track the number of genes used */
  int frame = 0;           /* Reading frame */
  int rbs_value = 0;       /* Best motif between exact and mismatch */
  int best_rbs[3] = {0};   /* RBS of the best gene in each frame */
  int best_type[3] = {0};  /* Type of the best gene in each frame */
  int best_index[3] = {0}; /* Index of best gene in each frame */
  double best_score[3] = {0};  /* Best scoring node in each frame */
  double wt = 0.0;         /* Shorthand for start weight to save space */
  double rbg[28] = {0};    /* RBS background - i.e. all nodes */
  double tbg[4] = {0};     /* Type background - i.e. all nodes */
  double dbg[16] = {0};    /* Dimer background - i.e. all nodes */
  double ubg[30][4] = {{0}};    /* Pair composition background */
  double rreal[28] = {0};  /* RBS foreground - i.e. only real genes */
  double treal[4] = {0};   /* Type foreground - i.e. only real genes */
  double dreal[16] = {0};  /* Dimer foreground - i.e. only real genes */
  double ureal[30][4] = {{0}};  /* Pair comp foreground - i.e. only real genes */
  double node_score = 0.0; /* Score for current node */
  double score_thresh = MIN_TRAIN_GENE_SCORE;
                           /* Minimum score for genes to be "real" */

  wt = train_data->start_weight;
  zero_start_weights(train_data, 0);

  /* Calculate the backgrounds that don't change */
  calc_start_background(seq, rseq, seq_length, train_data->trans_table, tbg);
  get_word_counts(2, seq, rseq, seq_length, dbg);
  normalize_array(dbg, 16);
  calc_upstream_background(nodes, num_nodes, ubg);

  /* Iterate SD_ITER times through the list of nodes                    */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few      */
  /* extra to be safe).                                                 */
  for (i = 0; i < SD_ITER; i++)
  {
    num_genes = 0;

    /* Recalculate the RBS motif background, since rbs weights change. */
    calc_sd_rbs_background(nodes, num_nodes, train_data->rbs_wt, rbg);

    /* Set real values to 0 since we're recounting */
    memset(rreal, 0, 28*sizeof(double));
    memset(treal, 0, 3*sizeof(double));

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
          num_genes++;
          rreal[best_rbs[frame]] += 1.0;
          treal[best_type[frame]] += 1.0;
          if (i == SD_ITER-1)
          {
            if (nodes[best_index[frame]].dimer >= 0)
            {
              dreal[nodes[best_index[frame]].dimer] += 1.0;
            }
            for (k = 0; k < 30; k++)
            {
              ureal[k][nodes[best_index[frame]].ups[k]] += 1.0;
            }
          }
        }
        best_score[frame] = 0.0;
        best_index[frame] = -1;
        best_rbs[frame] = 0;
        best_type[frame] = 0;
      }
      else if (nodes[j].strand == 1)
      {
        rbs_value = get_rbs_value(&nodes[j], train_data->rbs_wt);
        node_score = nodes[j].cscore + wt*train_data->rbs_wt[rbs_value] +
                     wt*train_data->type_wt[nodes[j].subtype];
        if (node_score >= best_score[frame])
        {
          best_score[frame] = node_score;
          best_index[frame] = j;
          best_type[frame] = nodes[j].subtype;
          best_rbs[frame] = rbs_value;
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
        if (best_score[frame] >= score_thresh &&
            nodes[best_index[frame]].index%3 == frame)
        {
          num_genes++;
          rreal[best_rbs[frame]] += 1.0;
          treal[best_type[frame]] += 1.0;
          if (i == SD_ITER-1)
          {
            if (nodes[j].dimer >= 0)
            {
              dreal[nodes[best_index[frame]].dimer] += 1.0;
            }
            for (k = 0; k < 30; k++)
            {
              ureal[k][nodes[best_index[frame]].ups[k]] += 1.0;
            }
          }
        }
        best_score[frame] = 0.0;
        best_index[frame] = -1;
        best_rbs[frame] = 0;
        best_type[frame] = 0;
      }
      else if (nodes[j].strand == -1)
      {
        rbs_value = get_rbs_value(&nodes[j], train_data->rbs_wt);
        node_score = nodes[j].cscore + wt*train_data->rbs_wt[rbs_value] +
                     wt*train_data->type_wt[nodes[j].subtype];
        if (node_score >= best_score[frame])
        {
          best_score[frame] = node_score;
          best_index[frame] = j;
          best_type[frame] = nodes[j].subtype;
          best_rbs[frame] = rbs_value;
        }
      }
    }

    /* Log score conversions */
    normalize_array(rreal, 28);
    create_log_score(train_data->rbs_wt, rreal, rbg, 28);
    normalize_array(treal, 4);
    create_log_score(train_data->type_wt, treal, tbg, 4);
    if (num_genes <= (double)num_nodes/SMALL_TRAIN_SET_NODES)
    {
      score_thresh /= 2.0;
    }
  }
  normalize_array(dreal, 16);
  create_log_score(train_data->dimer_wt, dreal, dbg, 16);
  for (i = 0; i < 30; i++)
  {
    normalize_array(ureal[i], 4);
    create_log_score(train_data->ups_wt[i], ureal[i], ubg[i], 4);
  }
/*
printf("base/real/bg/log\t%d\n", i);
for(j = 0; j < 4; j++) { printf("\t%d\t%.4f\t%.4f\t%.4f\n", j, ureal[i][j], ubg[i][j], train_data->ups_wt[i][j]); }
exit(0);
*/
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
  double zbg = 0.0;        /* Background for zero motif */
  double zreal = 0.0;      /* Zero motif weight for real genes */
  double ngenes = 0.0;     /* Gene counter, double in case of weighting */
  double score_thresh = MIN_TRAIN_GENE_SCORE;
                           /* Minimum score for genes to be "real" */

  wt = train_data->start_weight;
  zero_start_weights(train_data, 1);
  calc_type_background(nodes, num_nodes, tbg);

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

  for (i = 0; i < 4; i++)
  {
    train_data->type_wt[i] = 0.0;
    train_data->stop_wt[i] = 0.0;
  }
  for (i = 0; i < 16; i++)
  {
    train_data->dimer_wt[i] = 0.0;
  }
  for (i = 0; i < 30; i++)
  {
    for (j = 0; j < 4; j++)
    {
      train_data->ups_wt[i][j] = 0.0;
    }
  }

  if (which == 0)
  {
    for (i = 0; i < 28; i++)
    {
      train_data->rbs_wt[i] = 0.0;
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

  memset(tbg, 0, 4*sizeof(double));
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP || nodes[i].edge == 1)
    {
      continue;
    }
    tbg[nodes[i].subtype] += 1.0;
  }
  normalize_array(tbg, 4);
}

/* Calculate ATG/GTG/TTG background across whole sequence */
void calc_start_background(unsigned char *seq, unsigned char *rseq,
                           int seq_length, int tt, double *sbg)
{
  int i = 0;
  int type = 0;

  memset(sbg, 0, 4*sizeof(double));
  for (i = 0; i <= seq_length-3; i++)
  {
    type = get_start_type(seq, i, tt);
    if (type != -1)
    {
      sbg[type] += 1.0;
    }
    type = get_start_type(rseq, i, tt);
    if (type != -1)
    {
      sbg[type] += 1.0;
    }
  }
  normalize_array(sbg, 4);
}

/* Calculate TAA/TGA/TAG background across whole sequence */
void calc_stop_background(unsigned char *seq, unsigned char *rseq,
                          int seq_length, int tt, double *sbg)
{
  int i = 0;
  int type = 0;

  memset(sbg, 0, 4*sizeof(double));
  for (i = 0; i <= seq_length-3; i++)
  {
    type = get_stop_type(seq, i, tt);
    if (type != -1)
    {
      sbg[type] += 1.0;
    }
    type = get_stop_type(rseq, i, tt);
    if (type != -1)
    {
      sbg[type] += 1.0;
    }
  }
  normalize_array(sbg, 4);
}

/* Calculate the background for rbs given a set of weights */
void calc_sd_rbs_background(struct _node *nodes, int num_nodes, double *rbs_wt,
                            double *rbg)
{
  int i = 0;
  int rbs_switch = 0;      /* Best motif between exact and mismatch */

  memset(rbg, 0, 28*sizeof(double));
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP || nodes[i].edge == 1)
    {
      continue;
    }
    if (rbs_wt[nodes[i].rbs[0]] > rbs_wt[nodes[i].rbs[1]] + 1.0 ||
        nodes[i].rbs[1] == 0)
    {
      rbs_switch = nodes[i].rbs[0];
    }
    else if (rbs_wt[nodes[i].rbs[0]] < rbs_wt[nodes[i].rbs[1]] - 1.0 ||
             nodes[i].rbs[0] == 0)
    {
      rbs_switch = nodes[i].rbs[1];
    }
    else
    {
      rbs_switch = (int)dmax(nodes[i].rbs[0], nodes[i].rbs[1]);
    }
    rbg[rbs_switch] += 1.0;
  }
  normalize_array(rbg, 28);
}

/* Calculate the upstream background */
void calc_upstream_background(struct _node *nodes, int num_nodes, double ubg[30][4])
{
  int i = 0;
  int j = 0;

  /* Build the background of pair composition */
  for (i = 0; i < 30; i++)
  {
    for (j = 0; j < 4; j++)
    {
      ubg[i][j] = 0.0;
    }
  }

  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP || nodes[i].edge == 1)
    {
      continue;
    }
    for (j = 0; j < 30; j++)
    {
      ubg[j][nodes[i].ups[j]] += 1.0;
    }
  }
  for (i = 0; i < 30; i++)
  {
    normalize_array(ubg[i], 4);
  }
}

/* Make the components of a double array sum to 1.0 */
void normalize_array(double *arr, int size)
{
  int i = 0;
  double denom = 0.0;

  for (i = 0; i < size; i++)
  {
    denom += arr[i];
  }
  if (denom > 0.0)
  {
    for (i = 0; i < size; i++)
    {
      arr[i] /= denom;
    }
  }
}

/* Calculate a log-likelihood score for arr 1 and arr2 and store in arr3 */
void create_log_score(double *arr3, double *arr1, double *arr2, int size)
{
  int i = 0;
  double max = 4.0;
  double min = -4.0;

  if (arr1 == NULL || arr2 == NULL)
  {
    for (i = 0; i < size; i++)
    {
      if (arr3[i] == 0.0)
      {
        arr3[i] = min;
      }
      else
      {
        arr3[i] = log(size*arr3[i]);
      }
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      if (arr1[i] == 0.0 || arr2[i] == 0.0)
      {
        arr3[i] = min;
      }
      else
      {
        arr3[i] = log(arr1[i]/arr2[i]);
      }
    }
  }
  for (i = 0; i < size; i++)
  {
    if (arr3[i] > max)
    {
      arr3[i] = max;
    }
    if (arr3[i] < min)
    {
      arr3[i] = min;
    }
  }
}
