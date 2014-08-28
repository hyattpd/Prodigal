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

#include "node.h"

/******************************************************************************
  Adds nodes to the node list.  Genes must be >=90bp in length, unless they
  run off the edge, in which case they only have to be 60bp.
******************************************************************************/
int add_nodes(unsigned char *seq, unsigned char *rseq, unsigned char *useq,
              int seq_length, struct _node *nodes, int closed, int gap_mode,
              int trans_table)
{
  int i = 0;
  int num_nodes = 0;          /* Number of starts/stops */
  int last_stop[3] = {0};     /* Location of last stop seen in each frame */
  int saw_start[3] = {0};     /* 0/1 if we've seen a start in each frame */
  int min_dist[3] = {0};      /* Minimum distance required for gene */
                              /* Can be different for edge vs. normal */
  int edge[3] = {0};          /* 0/1 if the node in frame runs off edge */
  int sl_mod = 0;             /* Variable to help with frame math */
  int stop_codon = 0;         /* 0/1 if current codon is a stop codon */
  int start_codon = 0;        /* 0/1 if current codon is a start codon */

  /* Forward strand nodes */
  sl_mod = seq_length%3;
  for (i = 0; i < 3; i++)
  {
    last_stop[(i+sl_mod) % 3] = seq_length+i;
    while (last_stop[(i+sl_mod) % 3] + 2 > seq_length-1)
    {
      last_stop[(i+sl_mod) % 3] -= 3;
    }
  }
  /* Set edge stops */
  for (i = 0; i < 3; i++)
  {
    saw_start[i%3] = 0;
    if (is_stop(seq, last_stop[i%3], trans_table) == 1)
    {
      min_dist[i%3] = MIN_GENE;
      edge[i%3] = 0;
    }
    else
    {
      min_dist[i%3] = MIN_EDGE_GENE;
      edge[i%3] = 1;
    }
  }
  /* Work backwards through forward strand. */
  for (i = seq_length-3; i >= 0; i--)
  {
    stop_codon = is_stop(seq, i, trans_table);
    /* Stop Codon or a Run of N's to the Right */
    if (stop_codon == 1 || (i <= seq_length-12 && gap_mode < 2 &&
        gap_to_right(useq, i) == 1))
    {
      /* If previous stop had a valid start codon, we add it */
      if (saw_start[i%3] == 1)
      {
        nodes[num_nodes].edge = edge[i%3];
        nodes[num_nodes].index = last_stop[i%3];
        nodes[num_nodes].type = STOP;
        if (edge[i%3] == 1)
        {
          nodes[num_nodes].subtype = EDGE;
        }
        else
        {
          nodes[num_nodes].subtype = assign_stop_value(seq, last_stop[i%3]);
        }
        nodes[num_nodes].strand = 1;
        nodes[num_nodes++].stop_val = i;
      }
      /* Update minimum distances depending if edge node or not */
      if (stop_codon == 1)
      {
        min_dist[i%3] = MIN_GENE;
        edge[i%3] = 0;
      }
      else
      {
        min_dist[i%3] = MIN_EDGE_GENE;
        edge[i%3] = 1;
      }
      last_stop[i%3] = i;
      saw_start[i%3] = 0;
      continue;
    }
    /* No starts allowed if -c and edge of sequence or -z 1 and gap */
    if (edge[i%3] == 1 &&
        ((last_stop[i%3] + 5 > seq_length-1 && closed == 1) ||
         (last_stop[i%3] + 5 <= seq_length-1 && gap_mode == 1)))
    {
      continue;
    }

    /* Start Nodes */
    start_codon = is_start(seq, i, trans_table);
    /* Actual Start Codon */
    if (start_codon == 1 && ((last_stop[i%3]-i+3) >= min_dist[i%3]))
    {
      nodes[num_nodes].index = i;
      nodes[num_nodes].type = START;
      nodes[num_nodes].subtype = assign_start_value(seq, i);
      nodes[num_nodes].dimer = assign_dimer_value(seq, i);
      count_pair_composition(seq, seq_length, 1, i, nodes[num_nodes].pairs);
      normalize_array(nodes[num_nodes].pairs, 10);
      saw_start[i%3] = 1;
      nodes[num_nodes].stop_val = last_stop[i%3];
      nodes[num_nodes++].strand = 1;
    }
    /* Edge Start: Gap to the left with -z 0 or no -c and left edge */
    else if ((last_stop[i%3]-i) > MIN_EDGE_GENE &&
             ((i <= 2 && closed == 0) ||
             (gap_mode == 0 && i >= 9 && codon_has_n(useq, i) == 0 &&
             gap_to_left(useq, i) == 1)))
    {
      nodes[num_nodes].index = i;
      nodes[num_nodes].type = START;
      nodes[num_nodes].subtype = EDGE;
      saw_start[i%3] = 1;
      nodes[num_nodes].edge = 1;
      nodes[num_nodes].stop_val = last_stop[i%3];
      nodes[num_nodes++].strand = 1;
    }
  }
  /* Take care of final set of leftmost stop nodes */
  for (i = 0; i < 3; i++)
  {
    if (saw_start[i%3] == 1)
    {
      nodes[num_nodes].edge = edge[i%3];
      nodes[num_nodes].index = last_stop[i%3];
      nodes[num_nodes].type = STOP;
      if (edge[i%3] == 1)
      {
        nodes[num_nodes].subtype = EDGE;
      }
      else
      {
        nodes[num_nodes].subtype = assign_stop_value(seq, last_stop[i%3]);
      }
      nodes[num_nodes].strand = 1;
      nodes[num_nodes++].stop_val = i-6;
    }
  }

  /* Reverse strand nodes */
  for (i = 0; i < 3; i++)
  {
    last_stop[(i+sl_mod)%3] = seq_length+i;
    while (last_stop[(i+sl_mod)%3]+2 > seq_length-1)
    {
      last_stop[(i+sl_mod)%3]-=3;
    }
  }
  /* Set edge stops */
  for (i = 0; i < 3; i++)
  {
    saw_start[i%3] = 0;
    if (is_stop(rseq, last_stop[i%3], trans_table) == 1)
    {
      min_dist[i%3] = MIN_GENE;
      edge[i%3] = 1;
    }
    else
    {
      min_dist[i%3] = MIN_EDGE_GENE;
      edge[i%3] = 1;
    }
  }
  /* Work backwards through reverse strand. */
  for (i = seq_length-3; i >= 0; i--)
  {
    /* Stop Codon or a Run of N's to the Left (since useq is flipped) */
    stop_codon = is_stop(rseq, i, trans_table);
    if (stop_codon == 1 || (i <= seq_length-12 && gap_mode < 2 &&
        gap_to_left(useq, seq_length-3-i) == 1))
    {
      /* If previous stop had a valid start codon, we add it */
      if (saw_start[i%3] == 1)
      {
        nodes[num_nodes].edge = edge[i%3];
        nodes[num_nodes].index = seq_length-last_stop[i%3]-1;
        nodes[num_nodes].type = STOP;
        if (edge[i%3] == 1)
        {
          nodes[num_nodes].subtype = EDGE;
        }
        else
        {
          nodes[num_nodes].subtype = assign_stop_value(rseq, last_stop[i%3]);
        }
        nodes[num_nodes].strand = -1;
        nodes[num_nodes++].stop_val = seq_length-i-1;
      }
      /* Update minimum distances depending if edge node or not */
      if (is_stop(rseq, i, trans_table) == 1)
      {
        min_dist[i%3] = MIN_GENE;
        edge[i%3] = 0;
      }
      else
      {
        min_dist[i%3] = MIN_EDGE_GENE;
        edge[i%3] = 1;
      }
      last_stop[i%3] = i;
      saw_start[i%3] = 0;
      continue;
    }
    /* No starts allowed if -c and edge of sequence or -z 1 and gap */
    if (edge[i%3] == 1 &&
        ((last_stop[i%3] + 5 > seq_length-1 && closed == 1) ||
         (last_stop[i%3] + 5 <= seq_length-1 && gap_mode == 1)))
    {
      continue;
    }

    /* Start Nodes */
    start_codon = is_start(rseq, i, trans_table);
    /* Actual Start Codon */
    if (start_codon == 1 && ((last_stop[i%3]-i+3) >= min_dist[i%3]))
    {
      nodes[num_nodes].index = seq_length-i-1;
      nodes[num_nodes].type = START;
      nodes[num_nodes].subtype = assign_start_value(rseq, i);
      nodes[num_nodes].dimer = assign_dimer_value(rseq, i);
      count_pair_composition(rseq, seq_length, -1, seq_length-i-1,
                             nodes[num_nodes].pairs);
      normalize_array(nodes[num_nodes].pairs, 10);
      saw_start[i%3] = 1;
      nodes[num_nodes].stop_val = seq_length-last_stop[i%3]-1;
      nodes[num_nodes++].strand = -1;
    }
    /* Edge Start: Gap to the right (useq flipped) with -z 0 or no -c */
    /* and left edge */
    else if ((last_stop[i%3]-i) > MIN_EDGE_GENE &&
             ((i <= 2 && closed == 0) ||
             (gap_mode == 0 && i >= 9 &&
             codon_has_n(useq, seq_length-3-i) == 0 &&
             gap_to_right(useq, seq_length-3-i) == 1)))
    {
      nodes[num_nodes].index = seq_length-i-1;
      nodes[num_nodes].type = START;
      nodes[num_nodes].subtype = EDGE;
      saw_start[i%3] = 1;
      nodes[num_nodes].edge = 1;
      nodes[num_nodes].stop_val = seq_length-last_stop[i%3]-1;
      nodes[num_nodes++].strand = -1;
    }
  }
  /* Take care of final set of leftmost stop nodes */
  for (i = 0; i < 3; i++)
  {
    if (saw_start[i%3] == 1)
    {
      nodes[num_nodes].edge = edge[i%3];
      nodes[num_nodes].index = seq_length - last_stop[i%3] - 1;
      nodes[num_nodes].type = STOP;
      if (edge[i%3] == 1)
      {
        nodes[num_nodes].subtype = EDGE;
      }
      else
      {
        nodes[num_nodes].subtype = assign_stop_value(rseq, last_stop[i%3]);
      }
      nodes[num_nodes].strand = -1;
      nodes[num_nodes++].stop_val = seq_length-i+5;
    }
  }
  return num_nodes;
}

/* Memset nodes to 0 and return 0 */
void zero_nodes(struct _node *nodes, int num_nodes)
{
  memset(nodes, 0, num_nodes * sizeof(struct _node));
}

/* Check node allocation and realloc memory if we */
/* need more space. */
void check_node_allocation(struct _node **nodes, int seq_length)
{
  /* Grab more memory if sequence is larger than our default allocation */
  if (seq_length > STT_NOD*8)
  {
    *nodes = (struct _node *)
              realloc(*nodes, (int)(seq_length/8) * sizeof(struct _node));
    if (*nodes == NULL)
    {
      perror("Realloc failed on nodes\n\n");
      exit(11);
    }
  }
  zero_nodes(*nodes, seq_length/8);
}

/* Simple routine to zero out the node scores */
void reset_node_scores(struct _node *nodes, int num_nodes)
{
  int i = 0;
  int j = 0;
  for (i = 0; i < num_nodes; i++)
  {
    for (j = 0; j < 3; j++)
    {
      nodes[i].start_ptr[j] = 0;
      nodes[i].gc_score[j] = 0.0;
    }
    for (j = 0; j < 2; j++)
    {
      nodes[i].rbs[j] = 0;
    }
    nodes[i].score = 0.0;
    nodes[i].cscore = 0.0;
    nodes[i].sscore = 0.0;
    nodes[i].rscore = 0.0;
    nodes[i].tscore = 0.0;
    nodes[i].uscore = 0.0;
    nodes[i].trace_back = -1;
    nodes[i].trace_forward = -1;
    nodes[i].overlap_frame = -1;
    nodes[i].status = 0;
    nodes[i].gc_bias = 0;
    memset(&nodes[i].mot, 0, sizeof(struct _motif));
  }
}

/******************************************************************************
  Since dynamic programming can't go 'backwards', we have to record
  information about overlapping genes in order to build the models.  So, for
  example, in cases like 5'->3', 5'-3' overlapping on the same strand, we
  record information about the 2nd 5' end under the first 3' end's
  information.  For every stop, we calculate and store all the best starts
  that could be used in genes that overlap that 3' end.
******************************************************************************/
void record_overlapping_starts(struct _node *nodes, int num_nodes,
                               double start_weight, int stage)
{
  int i = 0;
  int j = 0;
  double max_score = 0.0;     /* Maximum score for an overlapping start */

  for (i = 0; i < num_nodes; i++)
  {
    for (j = 0; j < 3; j++)
    {
      nodes[i].start_ptr[j] = -1;
    }
    if (nodes[i].type == START || nodes[i].edge == 1)
    {
      continue;
    }
    if (nodes[i].strand == 1)
    {
      max_score = -100;
      for (j = i+3; j >= 0; j--)
      {
        if (j >= num_nodes || nodes[j].index > nodes[i].index+2)
        {
          continue;
        }
        if (nodes[j].index + MAX_SAM_OVLP < nodes[i].index)
        {
          break;
        }
        if (nodes[j].strand == 1 && nodes[j].type == START)
        {
          if (nodes[j].stop_val <= nodes[i].index)
          {
            continue;
          }
          if (stage == 0 && nodes[i].start_ptr[(nodes[j].index)%3] == -1)
          {
            nodes[i].start_ptr[(nodes[j].index)%3] = j;
          }
          else if (stage == 1 && (nodes[j].cscore + nodes[j].sscore +
                   intergenic_mod(&nodes[i], &nodes[j], start_weight) >
                   max_score))
          {
            nodes[i].start_ptr[(nodes[j].index)%3] = j;
            max_score = nodes[j].cscore + nodes[j].sscore +
                        intergenic_mod(&nodes[i], &nodes[j], start_weight);
          }
        }
      }
    }
    else
    {
      max_score = -100;
      for (j = i-3; j < num_nodes; j++)
      {
        if (j < 0 || nodes[j].index < nodes[i].index-2)
        {
          continue;
        }
        if (nodes[j].index - MAX_SAM_OVLP > nodes[i].index)
        {
          break;
        }
        if (nodes[j].strand == -1 && nodes[j].type == START)
        {
          if (nodes[j].stop_val >= nodes[i].index)
          {
            continue;
          }
          if (stage == 0 && nodes[i].start_ptr[(nodes[j].index)%3] == -1)
          {
            nodes[i].start_ptr[(nodes[j].index)%3] = j;
          }
          else if (stage == 1 && (nodes[j].cscore + nodes[j].sscore +
                   intergenic_mod(&nodes[j], &nodes[i], start_weight) >
                   max_score))
          {
            nodes[i].start_ptr[(nodes[j].index)%3] = j;
            max_score = nodes[j].cscore + nodes[j].sscore +
                        intergenic_mod(&nodes[j], &nodes[i], start_weight);
          }
        }
      }
    }
  }
}

/******************************************************************************
  This routine goes through all the ORFs and counts the relative frequency of
  the most common frame for G+C content.  In high GC genomes, this tends to be
  the third position.  In low GC genomes, this tends to be the first position.
  Genes will be selected as a training set based on the nature of this bias
  for this particular organism.  The "gc" variable contains the most common
  frame at a window centered on each position in the sequence.
******************************************************************************/
void frame_score(int *gc, struct _node *nodes, int num_nodes)
{
  int i = 0;
  int j = 0;
  int counter[3][3] = {{0}};     /* Counts for winning frame in each frame */
  int last[3] = {0};             /* Last stop codon we saw in each frame */
  int frame_mod = 0;             /* Relative-to-absolute frame modifier */
  int frame = 0;                 /* Frame of the current node */
  int best_frame = 0;            /* The final winning frame for a node */

  if (num_nodes == 0)
  {
    return;
  }
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      counter[i][j] = 0;
    }
  }
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
      counter[frame][(gc[nodes[i].index] + frame_mod)%3] = 1;
    }
    else if (nodes[i].strand == 1)
    {
      for (j = last[frame]-3; j >= nodes[i].index; j-=3)
      {
        counter[frame][(gc[j] + frame_mod)%3] ++;
      }
      best_frame = max_frame(counter[frame][0], counter[frame][1],
                             counter[frame][2]);
      nodes[i].gc_bias = best_frame;
      for (j = 0; j < 3; j++)
      {
        nodes[i].gc_score[j] = (3.0*counter[frame][j]);
        nodes[i].gc_score[j] /= (1.0*(nodes[i].stop_val - nodes[i].index + 3));
      }
      last[frame] = nodes[i].index;
    }
  }
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
      counter[frame][((3-gc[nodes[i].index]) + frame_mod)%3] = 1;
    }
    else if (nodes[i].strand == -1)
    {
      for (j = last[frame]+3; j <= nodes[i].index; j+=3)
      {
        counter[frame][((3-gc[j]) + frame_mod)%3]++;
      }
      best_frame = max_frame(counter[frame][0], counter[frame][1],
                             counter[frame][2]);
      nodes[i].gc_bias = best_frame;
      for (j = 0; j < 3; j++)
      {
        nodes[i].gc_score[j] = (3.0*counter[frame][j]);
        nodes[i].gc_score[j] /= (1.0*(nodes[i].index - nodes[i].stop_val + 3));
      }
      last[frame] = nodes[i].index;
    }
  }
}

/******************************************************************************
  Scoring function for all the start nodes.  This score has two factors:  (1)
  Coding, which is a composite of coding score and length, and (2) Start
  score, which is a composite of RBS score and ATG/TTG/GTG.
******************************************************************************/
void score_nodes(unsigned char *seq, unsigned char *rseq, int seq_length,
                 struct _node *nodes, int num_nodes,
                 struct _training *train_data, int closed, int mode)
{
  int i = 0;
  int j = 0;
  int num_edges = 0;            /* Number of open edges for partial gene */
  int min_anon_len = 0;         /* Min length for internal gene in anon mode */
  double neg_fac = 0.0;         /* Negative factor to apply to pos. scores */
  double pos_fac = 0.0;         /* Positive factor to apply to neg. scores */
  double rbs_exact = 0.0;       /* Best score for an exact SD motif */
  double rbs_mismatch = 0.0;    /* Best score for a SD motif w/mismatch */
  double sd_score = 0.0;        /* Max score of exact/mismatch rbs scores */

  /* Step 1: Calculate raw coding potential for every start-stop pair. */
  calc_orf_gc(seq, nodes, num_nodes);
  calc_coding_score(seq, rseq, seq_length, nodes, num_nodes,
                    train_data->trans_table, train_data->gc,
                    train_data->gene_dc);

  /* Step 2: Calculate raw RBS Scores for every start node. */
  if (train_data->uses_sd == 1)
  {
    sd_rbs_score(seq, rseq, seq_length, nodes, num_nodes, train_data->rbs_wt);
  }
  else
  {
    for (i = 0; i < num_nodes; i++)
    {
      if (nodes[i].type == STOP || nodes[i].edge == 1)
      {
        continue;
      }
      find_best_nonsd_motif(train_data, seq, rseq, seq_length, &nodes[i], 2);
    }
  }

  /* Step 3: Score the start nodes */
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP)
    {
      continue;
    }

    /* Does this gene run off the edge? */
    num_edges = 0;
    if (nodes[i].edge == 1)
    {
      num_edges++;
    }
    if ((nodes[i].strand == 1 && is_stop(seq, nodes[i].stop_val,
        train_data->trans_table) == 0) || (nodes[i].strand == -1 &&
        is_stop(rseq, seq_length-1-nodes[i].stop_val,
                train_data->trans_table) == 0))
    {
      num_edges++;
    }

    /* Edge Nodes : stops with no starts, give a small bonus */
    if (nodes[i].edge == 1)
    {
      nodes[i].tscore = EDGE_BONUS*train_data->start_weight/(double)num_edges;
      nodes[i].uscore = 0.0;
      nodes[i].rscore = 0.0;
    }

    else
    {

      /* Start Sequence Score: Upstream Dimer plus Codon */
      nodes[i].tscore = train_data->type_wt[nodes[i].subtype] *
                        train_data->start_weight;
      nodes[i].tscore += 2 * train_data->dimer_wt[nodes[i].dimer];

      /* RBS Motif Score */
      rbs_exact = train_data->rbs_wt[nodes[i].rbs[0]];
      rbs_mismatch = train_data->rbs_wt[nodes[i].rbs[1]];
      sd_score = dmax(rbs_exact, rbs_mismatch) * train_data->start_weight;
      if (train_data->uses_sd == 1)
      {
        nodes[i].rscore = sd_score;
      }
      else
      {
        nodes[i].rscore = train_data->start_weight*nodes[i].mot.score;
        if (nodes[i].rscore < sd_score && train_data->no_mot > -0.5)
        {
          nodes[i].rscore = sd_score;
        }
      }

      /* Secondary Structure Score */
      for (j = 0; j < 10; j++)
      {
        nodes[i].uscore += 4 * (SSTRUCT_SIZE - 6) * train_data->pair_wt[j] * nodes[i].pairs[j];
/* printf("node %d, pair %d, adding %.4f pairs with %.4f weight, new score %.4f\n", i, j, nodes[i].pairs[j], train_data->pair_wt[j], nodes[i].uscore); */
      }

      /****************************************************************
      ** Penalize upstream score if choosing this start would stop   **
      ** the gene from running off the edge.                         **
      ****************************************************************/
      if (closed == 0 && nodes[i].index <= 2 && nodes[i].strand == 1)
      {
        nodes[i].uscore += EDGE_UPS*train_data->start_weight;
      }
      else if (closed == 0 && nodes[i].index >= seq_length-3 &&
               nodes[i].strand == -1)
      {
        nodes[i].uscore += EDGE_UPS*train_data->start_weight;
      }
      else if (i < 500 && nodes[i].strand == 1)
      {
        for (j = i-1; j >= 0; j--)
        {
          if (nodes[j].edge == 1 && nodes[i].stop_val == nodes[j].stop_val)
          {
            nodes[i].uscore += EDGE_UPS*train_data->start_weight;
            break;
          }
        }
      }
      else if (i >= num_nodes-500 && nodes[i].strand == -1)
      {
        for (j = i+1; j < num_nodes; j++)
        {
          if (nodes[j].edge == 1 && nodes[i].stop_val == nodes[j].stop_val)
          {
            nodes[i].uscore += EDGE_UPS*train_data->start_weight;
            break;
          }
        }
      }
    }

    /* Convert starts at base 1 and seq_length to edge genes if closed = 0 */
    if (((nodes[i].index <= 2 && nodes[i].strand == 1) ||
        (nodes[i].index >= seq_length-3 && nodes[i].strand == -1)) &&
        nodes[i].edge == 0 && closed == 0)
    {
      num_edges++;
      nodes[i].edge = 1;
      nodes[i].tscore = 0.0;
      nodes[i].uscore = EDGE_BONUS*train_data->start_weight/(double)num_edges;
      nodes[i].rscore = 0.0;
    }

    /* Penalize starts with no stop codon */
    if (nodes[i].edge == 0 && num_edges == 1)
    {
      nodes[i].uscore -= 0.5*EDGE_BONUS*train_data->start_weight;
    }

    /* Penalize non-edge genes < 250bp */
    if (num_edges == 0 && abs(nodes[i].index-nodes[i].stop_val) < 250)
    {
      neg_fac = 250.0/(float)abs(nodes[i].index-nodes[i].stop_val);
      pos_fac = (float)abs(nodes[i].index-nodes[i].stop_val)/250.0;
      if (nodes[i].rscore < 0)
      {
        nodes[i].rscore *= neg_fac;
      }
      if (nodes[i].uscore < 0)
      {
        nodes[i].uscore *= neg_fac;
      }
      if (nodes[i].tscore < 0)
      {
        nodes[i].tscore *= neg_fac;
      }
      if (nodes[i].rscore > 0)
      {
        nodes[i].rscore *= pos_fac;
      }
      if (nodes[i].uscore > 0)
      {
        nodes[i].uscore *= pos_fac;
      }
      if (nodes[i].tscore > 0)
      {
        nodes[i].tscore *= pos_fac;
      }
    }

    /**************************************************************/
    /* Coding Penalization in Anonymous Fragments:    Internal    */
    /* genes must have a score of 5.0 and be >= 120bp.  High GC   */
    /* genes are also penalized.                                  */
    /**************************************************************/
    if (mode == MODE_ANON && seq_length < 3000 && num_edges == 0 &&
        (nodes[i].cscore < 5.0 || abs(nodes[i].index-nodes[i].stop_val < 120)))
    {
      nodes[i].cscore -= META_PEN*dmax(0, (3000-seq_length)/2700.0);
    }

    /* Base Start Score */
    nodes[i].sscore = nodes[i].tscore + nodes[i].rscore + nodes[i].uscore;

    /**************************************************************/
    /* Penalize starts if coding is negative.  Larger penalty for */
    /* edge genes, since the start is offset by a smaller amount  */
    /* of coding than normal.                                     */
    /**************************************************************/
    if (nodes[i].cscore < 0.0)
    {
      if (num_edges > 0 && nodes[i].edge == 0)
      {
        if (mode != MODE_ANON || seq_length > 1500)
        {
          nodes[i].sscore -= train_data->start_weight;
        }
        else
        {
          nodes[i].sscore -= (10.31 - 0.004*seq_length);
        }
      }
      else if (mode == MODE_ANON && seq_length < 3000 && nodes[i].edge == 1)
      {
        min_anon_len = sqrt(seq_length)*5.0;
        if (abs(nodes[i].index-nodes[i].stop_val) >= min_anon_len)
        {
          if (nodes[i].cscore >= 0)
          {
            nodes[i].cscore = -1.0;
          }
          nodes[i].sscore = 0.0;
          nodes[i].uscore = 0.0;
        }
      }
      else
      {
        nodes[i].sscore -= 0.5;
      }
    }
    else if (nodes[i].cscore < 5.0 && mode == MODE_ANON &&
             abs(nodes[i].index-nodes[i].stop_val < 120) &&
             nodes[i].sscore < 0.0)
    {
      nodes[i].sscore -= train_data->start_weight;
    }
  }
}

/* Calculate the GC Content for each start-stop pair */
void calc_orf_gc(unsigned char *seq, struct _node *nodes, int num_nodes)
{
  int i = 0;
  int j = 0;
  int last[3] = {0};         /* Last stop seen in this frame */
  int frame = 0;             /* Frame of current node position */
  double gc[3] = {0};        /* Current GC count in each frame */
  double gene_size = 0.0;    /* Current gene size to divide GC by */

  /* Go through each start-stop pair and calculate the %GC of the gene */
  for (i = 0; i < 3; i++)
  {
    gc[i] = 0.0;
  }
  for (i = num_nodes-1; i >= 0; i--)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == 1 && nodes[i].type == STOP)
    {
      last[frame] = nodes[i].index;
      gc[frame] = is_gc(seq, nodes[i].index) + is_gc(seq, nodes[i].index+1) +
                  is_gc(seq, nodes[i].index+2);
    }
    else if (nodes[i].strand == 1)
    {
      for (j = last[frame]-3; j >= nodes[i].index; j-=3)
      {
        gc[frame] += is_gc(seq, j) + is_gc(seq, j+1) + is_gc(seq, j+2);
      }
      gene_size = (float)(abs(nodes[i].stop_val-nodes[i].index)+3.0);
      nodes[i].gc_cont = gc[frame]/gene_size;
      last[frame] = nodes[i].index;
    }
  }
  for (i = 0; i < 3; i++)
  {
    gc[i] = 0.0;
  }
  for (i = 0; i < num_nodes; i++)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == -1 && nodes[i].type == STOP)
    {
      last[frame] = nodes[i].index;
      gc[frame] = is_gc(seq, nodes[i].index) + is_gc(seq, nodes[i].index-1) +
                  is_gc(seq, nodes[i].index-2);
    }
    else if (nodes[i].strand == -1)
    {
      for (j = last[frame]+3; j <= nodes[i].index; j+=3)
      {
        gc[frame] += is_gc(seq, j) + is_gc(seq, j+1) + is_gc(seq, j+2);
      }
      gene_size = (float)(abs(nodes[i].stop_val-nodes[i].index)+3.0);
      nodes[i].gc_cont = gc[frame]/gene_size;
      last[frame] = nodes[i].index;
    }
  }
}

/******************************************************************************
  Score each candidate's coding.  We also sharpen coding/noncoding thresholds
  to prevent choosing interior starts when there is strong coding continuing
  upstream.
******************************************************************************/
void calc_coding_score(unsigned char *seq, unsigned char *rseq, int seq_length,
                       struct _node *nodes, int num_nodes, int trans_table,
                       double gc, double *hex_probs)
{
  int i = 0;
  int j = 0;
  int last[3] = {0};           /* Last index seen in this frameame */
  int frame = 0;               /* Frame of current node */
  double score[3] = {0};       /* Running coding score in each frameame */
  double len_fac = 0.0;        /* Length factor added for long genes */
  double prob_no_stop = 0.0;   /* Prob of seeing a non-stop-codon */
  double gene_size = 0.0;      /* Size of current gene */

  prob_no_stop = 1.0-prob_stop(trans_table, gc);

  /* Initial Pass: Score coding potential (start->stop) */
  for (i = 0; i < 3; i++)
  {
    score[i] = 0.0;
  }
  for (i = num_nodes-1; i >= 0; i--)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == 1 && nodes[i].type == STOP)
    {
      last[frame] = nodes[i].index - COD_SKIP;
      score[frame] = 0.0;
    }
    else if (nodes[i].strand == 1)
    {
      for (j = last[frame]-3; j >= nodes[i].index; j-=3)
      {
        score[frame] += hex_probs[mer_index(6, seq, j + COD_SKIP)];
      }
      nodes[i].cscore = score[frame];
      last[frame] = nodes[i].index;
    }
  }
  for (i = 0; i < 3; i++)
  {
    score[i] = 0.0;
  }
  for (i = 0; i < num_nodes; i++)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == -1 && nodes[i].type == STOP)
    {
      last[frame] = nodes[i].index + COD_SKIP;
      score[frame] = 0.0;
    }
    else if (nodes[i].strand == -1)
    {
      for (j = last[frame]+3; j <= nodes[i].index; j+=3)
      {
        score[frame] +=
          hex_probs[mer_index(6, rseq, seq_length - j - 1 + COD_SKIP)];
      }
      nodes[i].cscore = score[frame];
      last[frame] = nodes[i].index;
    }
  }

  /* Second Pass: Penalize start nodes with ascending coding upstream */
  for (i = 0; i < 3; i++)
  {
    score[i] = -10000.0;
  }
  for (i = 0; i < num_nodes; i++)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == 1 && nodes[i].type == STOP)
    {
      score[frame] = -10000.0;
    }
    else if (nodes[i].strand == 1)
    {
      if (nodes[i].cscore > score[frame])
      {
        score[frame] = nodes[i].cscore;
      }
      else
      {
        nodes[i].cscore -= (score[frame] - nodes[i].cscore);
      }
    }
  }
  for (i = 0; i < 3; i++)
  {
    score[i] = -10000.0;
  }
  for (i = num_nodes-1; i >= 0; i--)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == -1 && nodes[i].type == STOP)
    {
      score[frame] = -10000.0;
    }
    else if (nodes[i].strand == -1)
    {
      if (nodes[i].cscore > score[frame])
      {
        score[frame] = nodes[i].cscore;
      }
      else
      {
        nodes[i].cscore -= (score[frame] - nodes[i].cscore);
      }
    }
  }

  /* Third Pass: Add length-based factor to the score      */
  /* Penalize start nodes based on length to their left    */
  for (i = 0; i < num_nodes; i++)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == 1 && nodes[i].type == STOP)
    {
      score[frame] = -10000.0;
    }
    else if (nodes[i].strand == 1)
    {
      gene_size = ((float)(abs(nodes[i].stop_val-nodes[i].index) + 3.0))/3.0;
      if (gene_size > 1000.0)
      {
        len_fac = log((1-pow(prob_no_stop, 1000.0))/pow(prob_no_stop, 1000.0));
        len_fac -= log((1-pow(prob_no_stop, 80))/pow(prob_no_stop, 80));
        len_fac *= (gene_size - 80) / 920.0;
      }
      else
      {
        len_fac = log((1-pow(prob_no_stop, gene_size))/
                      pow(prob_no_stop, gene_size));
        len_fac -= log((1-pow(prob_no_stop, 80))/pow(prob_no_stop, 80));
      }
      if (len_fac > score[frame])
      {
        score[frame] = len_fac;
      }
      else
      {
        len_fac -= dmax(dmin(score[frame] - len_fac, len_fac), 0);
      }
      if (len_fac > 3.0 && nodes[i].cscore < 0.5*len_fac)
      {
        nodes[i].cscore = 0.5*len_fac;
      }
      nodes[i].cscore += len_fac;
    }
  }
  for (i = num_nodes-1; i >= 0; i--)
  {
    frame = (nodes[i].index)%3;
    if (nodes[i].strand == -1 && nodes[i].type == STOP)
    {
      score[frame] = -10000.0;
    }
    else if (nodes[i].strand == -1)
    {
      gene_size = ((float)(abs(nodes[i].stop_val-nodes[i].index)+3.0))/3.0;
      if (gene_size > 1000.0)
      {
        len_fac = log((1-pow(prob_no_stop, 1000.0))/pow(prob_no_stop, 1000.0));
        len_fac -= log((1-pow(prob_no_stop, 80))/pow(prob_no_stop, 80));
        len_fac *= (gene_size - 80) / 920.0;
      }
      else
      {
        len_fac = log((1-pow(prob_no_stop, gene_size))/
                      pow(prob_no_stop, gene_size));
        len_fac -= log((1-pow(prob_no_stop, 80))/pow(prob_no_stop, 80));
      }
      if (len_fac > score[frame])
      {
        score[frame] = len_fac;
      }
      else
      {
        len_fac -= dmax(dmin(score[frame] - len_fac, len_fac), 0);
      }
      if (len_fac > 3.0 && nodes[i].cscore < 0.5*len_fac)
      {
        nodes[i].cscore = 0.5*len_fac;
      }
      nodes[i].cscore += len_fac;
    }
  }
}

/******************************************************************************
  Shine-Dalgarno RBS Scoring Function: Calculate the best SD motif and then
  multiply it by the appropriate weight for that motif (determined in the
  start training function).
******************************************************************************/
void sd_rbs_score(unsigned char *seq, unsigned char *rseq, int seq_length,
                  struct _node *nodes, int num_nodes, double *rbs_weight)
{
  int i = 0;
  int j = 0;
  int motif[2] = {0};      /* Bins of best exact and mismatch RBS motifs */

  /* Scan all starts looking for RBS's */
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP || nodes[i].edge == 1)
    {
      continue;
    }
    nodes[i].rbs[0] = 0;
    nodes[i].rbs[1] = 0;
    if (nodes[i].strand == 1)
    {
      for (j = nodes[i].index - 20; j <= nodes[i].index - 6; j++)
      {
        if (j < 0)
        {
          continue;
        }
        motif[0] = shine_dalgarno_exact(seq, j, nodes[i].index, rbs_weight);
        motif[1] = shine_dalgarno_mismatch(seq, j, nodes[i].index, rbs_weight);
        if (motif[0] > nodes[i].rbs[0])
        {
          nodes[i].rbs[0] = motif[0];
        }
        if (motif[1] > nodes[i].rbs[1])
        {
          nodes[i].rbs[1] = motif[1];
        }
      }
    }
    else if (nodes[i].strand == -1)
    {
      for (j = seq_length-nodes[i].index-21; j <= seq_length-nodes[i].index-7;
           j++)
      {
        if (j > seq_length-1)
        {
          continue;
        }
        motif[0] = shine_dalgarno_exact(rseq, j, seq_length-1-nodes[i].index,
                                        rbs_weight);
        motif[1] = shine_dalgarno_mismatch(rseq, j, seq_length - 1 -
                                           nodes[i].index, rbs_weight);
        if (motif[0] > nodes[i].rbs[0])
        {
          nodes[i].rbs[0] = motif[0];
        }
        if (motif[1] > nodes[i].rbs[1])
        {
          nodes[i].rbs[1] = motif[1];
        }
      }
    }
  }
}

/******************************************************************************
  Given the weights for various motifs/distances from the training file,
  return the highest scoring mer/spacer combination of 3-6bp motifs with a
  spacer ranging from 3bp to 15bp.  In the final stage of start training, only
  good scoring motifs are returned.
******************************************************************************/
void find_best_nonsd_motif(struct _training *train_data, unsigned char *seq,
                           unsigned char *rseq, int seq_length,
                           struct _node *nodes, int stage)
{
  int i = 0;
  int j = 0;
  int start = 0;               /* Position of node in working sequence */
  int spacer = 0;              /* Spacer length */
  int spacer_index = 0;        /* Spacer index in weight array */
  int motif_index = 0;         /* Index of the motif (0-4095) */
  int max_spacer = 0;          /* Best motif spacer */
  int max_spacer_index = 0;    /* Best motif spacer index */
  int max_len = 0;             /* Best motif size */
  int max_index = 0;           /* Best motif index */
  double max_score = -100.0;   /* Best motif score */
  double score = 0.0;          /* Current motif score */
  unsigned char *wseq = NULL;  /* Working sequence - seq or rseq */

  if (nodes->type == STOP || nodes->edge == 1)
  {
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

  /* Work backwards from 6 base motifs down to 3 */
  /* If a motif has the best score, that's the one we go with */
  for (i = 3; i >= 0; i--)
  {
    for (j = start-18-i; j <= start-6-i; j++)
    {
      if (j < 0)
      {
        continue;
      }
      spacer = start-j-i-3;
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
      motif_index = mer_index(i+3, wseq, j);
      score = train_data->mot_wt[i][spacer_index][motif_index];
      if (score > max_score)
      {
        max_score = score;
        max_spacer_index = spacer_index;
        max_spacer = spacer;
        max_index = motif_index;
        max_len = i+3;
      }
    }
  }

  /* If in final stage, we don't accept bad scores (the score must be */
  /* at least a factor of 2 (0.69 log) better than not having any motif. */
  if (stage == 2 && (max_score == -4.0 || max_score < train_data->no_mot+0.69))
  {
    nodes->mot.index = 0;
    nodes->mot.len = 0;
    nodes->mot.spacer_index = 0;
    nodes->mot.spacer = 0;
    nodes->mot.score = train_data->no_mot;
  }
  else
  {
    nodes->mot.index = max_index;
    nodes->mot.len = max_len;
    nodes->mot.spacer_index = max_spacer_index;
    nodes->mot.spacer = max_spacer;
    nodes->mot.score = max_score;
  }
}

/* Return rbs value for a node */
int get_rbs_value(struct _node *n1, double *rbs_wt)
{
  if (rbs_wt[n1->rbs[0]] > rbs_wt[n1->rbs[1]] + 1.0 ||
      n1->rbs[1] == 0)
  {
    return n1->rbs[0];
  }
  else if (rbs_wt[n1->rbs[0]] < rbs_wt[n1->rbs[1]] - 1.0 ||
           n1->rbs[0] == 0)
  {
    return n1->rbs[1];
  }
  else
  {
    return (int)dmax(n1->rbs[0], n1->rbs[1]);
  }
}

/******************************************************************************
  When connecting two genes, we add a bonus for the -1 and -4 base overlaps on
  the same strand, which often signify an operon and negate the need for an
  RBS for the second gene.  In addition, we add a slight bonus when genes are
  close and a slight penalty when switching strands or having a large
  space between genes.
******************************************************************************/
double intergenic_mod(struct _node *n1, struct _node *n2, double start_weight)
{
  int dist = 0;
  int overlap = 0;         /* 1 = the nodes overlap, 0 = they don't */
  double ret_val = 0.0;    /* Return value for this function */

  if ((n1->strand == 1 && n2->strand == 1 &&
      (n1->index + 2 == n2->index || n1->index - 1 == n2->index)) ||
      (n1->strand == -1 && n2->strand == -1 &&
      (n1->index + 2 == n2->index || n1->index - 1 == n2->index)))
  {
    if (n1->strand == 1 && n2->rscore < 0)
    {
      ret_val -= n2->rscore;
    }
    if (n1->strand == -1 && n1->rscore < 0)
    {
      ret_val -= n1->rscore;
    }
    if (n1->strand == 1 && n2->uscore < 0)
    {
      ret_val -= n2->uscore;
    }
    if (n1->strand == -1 && n1->uscore < 0)
    {
      ret_val -= n1->uscore;
    }
  }
  dist = abs(n1->index-n2->index);
  if (n1->strand == 1 && n2->strand == 1 && n1->index+2 >= n2->index)
  {
    overlap = 1;
  }
  else if (n1->strand == -1 && n2->strand == -1 && n1->index >= n2->index+2)
  {
    overlap = 1;
  }
  if (dist > 3*OPER_DIST || n1->strand != n2->strand)
  {
    ret_val -= 0.15 * start_weight;
  }
  else if ((dist <= OPER_DIST && overlap == 0) || dist < 0.25*OPER_DIST)
  {
    ret_val += (2.0 - (double)(dist)/OPER_DIST) * 0.15 * start_weight;
  }
  return ret_val;
}

/* Default sorting routine for nodes */
int compare_nodes(const void *v1, const void *v2)
{
  struct _node *n1 = (struct _node *)v1;
  struct _node *n2 = (struct _node *)v2;
  if (n1->index < n2->index)
  {
    return -1;
  }
  if (n1->index > n2->index)
  {
    return 1;
  }
  if (n1->strand > n2->strand)
  {
    return -1;
  }
  if (n1->strand < n2->strand)
  {
    return 1;
  }
  return 0;
}

/* Sorts all nodes by common stop */
int stopcmp_nodes(const void *v1, const void *v2)
{
  struct _node *n1 = (struct _node *)v1;
  struct _node *n2 = (struct _node *)v2;
  if (n1->stop_val < n2->stop_val)
  {
    return -1;
  }
  if (n1->stop_val > n2->stop_val)
  {
    return 1;
  }
  if (n1->strand > n2->strand)
  {
    return -1;
  }
  if (n1->strand < n2->strand)
  {
    return 1;
  }
  if (n1->index < n2->index)
  {
    return -1;
  }
  if (n1->index > n2->index)
  {
    return 1;
  }
  return 0;
}
