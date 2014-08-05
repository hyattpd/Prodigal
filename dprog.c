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

#include "dprog.h"

/******************************************************************************
  Basic dynamic programming routine to find the optimal path of potential
  genes through the genome.  The routine is complicated by the fact genes
  can and frequently do overlap, but we've recorded information about
  optimal overlaps before executing this function.

  The 'stage' variable determines which of two scoring functions is used
  for the dynamic programming.
    stage = 0: The GC frame bias is used as the coding score.
    stage = 1: The hexamer coding scores are used instead.
******************************************************************************/
int dynamic_programming(struct _node *nodes, int num_nodes, double *frame_bias,
                        double start_weight, int stage)
{
  int i = 0;
  int j = 0;
  int low_bound = -1;            /* Index of lowest node to connect */
  int max_index = -1;            /* Pointer to highest scoring node */
  double max_score = -1.0;       /* Score of highest scoring node */

  /* Return -1 if there are no starts/stops in this sequence */
  if (num_nodes == 0)
  {
    return -1;
  }
  /* Initialize dynamic programming scores and trace pointers */
  for (i = 0; i < num_nodes; i++)
  {
    nodes[i].score = 0;
    nodes[i].trace_back = -1;
    nodes[i].trace_forward = -1;
  }
  /* Dynamic programming */
  for (i = 0; i < num_nodes; i++)
  {
    /* Get farthest node away from which to start connecting */
    low_bound = find_farthest_allowable_node(nodes, i);
    /* Score all connections from nodes between low_bound and the current
       node-1 with the current_node */
    for (j = low_bound; j < i; j++)
    {
      score_connection(nodes, j, i, frame_bias, start_weight, stage);
    }
  }
  /* Locate the highest scoring node in the traceback and record its index */
  for (i = num_nodes-1; i >= 0; i--)
  {
    if (nodes[i].strand == 1 && nodes[i].type == START)
    {
      continue;
    }
    if (nodes[i].strand == -1 && nodes[i].type == STOP)
    {
      continue;
    }
    if (nodes[i].score > max_score)
    {
      max_score = nodes[i].score;
      max_index = i;
    }
  }
  if (nodes[max_index].trace_back == -1)
  {
    return -1;
  }

  /* Insert nodes into the model missing due to overlaps. */
  untangle_overlaps(nodes, max_index);

  /* Dynamic programming did a traceback, but now we want to record
     trace-forward pointers, so we can get the genes in forward order. */
  record_trace_forward_pointers(nodes, max_index);

  /* Mark all nodes in the model as status 1 */
  mark_valid_gene_nodes(nodes, max_index);

  return max_index;
}

/******************************************************************************
  Dynamic programming is n log n, so we set up a distance constraint so that
  we will only try to connect two nodes within MAX_NODE_DIST of each other in
  the sorted list of nodes.  We make an exception for same-gene connections
  that represent giant ORFs.  This function returns the maximally distant
  node to start making connections to the current node from.
******************************************************************************/
int find_farthest_allowable_node(struct _node *nodes, int index)
{
  int low_bound = 0;

  if (index >= MAX_NODE_DIST)
  {
    low_bound = index-MAX_NODE_DIST;
  }
  if (nodes[index].strand == -1 && nodes[index].type == START &&
      nodes[low_bound].index >= nodes[index].stop_val)
  {
    while (low_bound >= 0 &&
           nodes[low_bound].index != nodes[index].stop_val)
    {
      low_bound--;
    }
  }
  if (nodes[index].strand == 1 && nodes[index].type == STOP &&
      nodes[low_bound].index >= nodes[index].stop_val)
  {
    while (low_bound >= 0 &&
           nodes[low_bound].index != nodes[index].stop_val)
    {
      low_bound--;
    }
  }
  if (low_bound < MAX_NODE_DIST)
  {
    low_bound = 0;
  }
  else
  {
    low_bound = low_bound-MAX_NODE_DIST;
  }
  return low_bound;
}

/******************************************************************************
  This routine scores the connection between two nodes, the most basic of which
  is 5'fwd->3'fwd (gene) and 3'rev->5'rev (rev gene).  If the connection ending
  at n2 is the maximal scoring model, it updates the pointers in the dynamic
  programming model.  n3 is used to handle overlaps, i.e. cases where 5->3'
  overlaps 5'->3' on the same strand. In this case, 3' connects directly to 3',
  and n3 is used to untangle the 5' end of the second gene.
******************************************************************************/
void score_connection(struct _node *nodes, int node_a, int node_b,
                      double *frame_bias, double start_weight, int stage)
{
  struct _node *n1 = &(nodes[node_a]);   /* Short reference to left node */
  struct _node *n2 = &(nodes[node_b]);   /* Short reference to right node */
  struct _node *n3 = NULL;               /* Reference to node inside overlap */
  int i = 0;
  int left = n1->index;                  /* Indices in nodes list */
  int right = n2->index;
  int bound = -1;                        /* Overlap-related tracking vars */
  int overlap = 0;
  int max_frame = -1;
  double max_val = -1000.0;
  double score = 0.0;                    /* Score of this connection */

  /***********************/
  /* Invalid Connections */
  /***********************/

  /* 5'fwd->5'fwd, 5'rev->5'rev */
  if (n1->type == START && n2->type == START &&
      n1->strand == n2->strand)
  {
    return;
  }

  /* 5'fwd->5'rev, 5'fwd->3'rev */
  else if (n1->strand == 1 && n1->type == START && n2->strand == -1)
  {
    return;
  }

  /* 3'rev->5'fwd, 3'rev->3'fwd) */
  else if (n1->strand == -1 && n1->type == STOP && n2->strand == 1)
  {
    return;
  }

  /* 5'rev->3'fwd */
  else if (n1->strand == -1 && n1->type == START && n2->strand == 1 &&
           n2->type == STOP)
  {
    return;
  }

  /******************/
  /* Edge Artifacts */
  /******************/
  if (n1->trace_back == -1 && n1->strand == 1 && n1->type == STOP)
  {
    return;
  }
  if (n1->trace_back == -1 && n1->strand == -1 && n1->type == START)
  {
    return;
  }

  /*********/
  /* Genes */
  /*********/

  /* 5'fwd->3'fwd */
  else if (n1->strand == n2->strand && n1->strand == 1 &&
           n1->type == START && n2->type == STOP)
  {
    if (n2->stop_val >= n1->index)
    {
      return;
    }
    if (n1->index % 3 != n2->index % 3)
    {
      return;
    }
    right += 2;
    if (stage == 0)
    {
      score = frame_bias[0]*n1->gc_score[0] + frame_bias[1]*n1->gc_score[1] +
              frame_bias[2]*n1->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n1->cscore + n1->sscore;
    }
  }

  /* 3'rev->5'rev */
  else if (n1->strand == n2->strand && n1->strand == -1 &&
           n1->type == STOP && n2->type == START)
  {
    if (n1->stop_val <= n2->index)
    {
      return;
    }
    if (n1->index % 3 != n2->index % 3)
    {
      return;
    }
    left -= 2;
    if (stage == 0)
    {
      score = frame_bias[0]*n2->gc_score[0] + frame_bias[1]*n2->gc_score[1] +
              frame_bias[2]*n2->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n2->cscore + n2->sscore;
    }
  }

  /********************************/
  /* Intergenic Space (Noncoding) */
  /********************************/

  /* 3'fwd->5'fwd */
  else if (n1->strand == 1 && n1->type == STOP && n2->strand == 1 &&
           n2->type == START)
  {
    left += 2;
    if (left >= right)
    {
      return;
    }
    if (stage == 1)
    {
      score = intergenic_mod(n1, n2, start_weight);
    }
  }

  /* 3'fwd->3'rev */
  else if (n1->strand == 1 && n1->type == STOP && n2->strand == -1 &&
           n2->type == STOP)
  {
    left += 2;
    right -= 2;
    if (left >= right)
    {
      return;
    }
    /* Overlapping Gene Case 2: Three consecutive overlapping genes f r r */
    max_frame = -1;
    max_val = 0.0;
    for (i = 0; i < 3; i++)
    {
      if (n2->star_ptr[i] == -1)
      {
        continue;
      }
      n3 = &(nodes[n2->star_ptr[i]]);
      overlap = left - n3->stop_val + 3;
      if (overlap <= 0 || overlap >= MAX_OPP_OVLP)
      {
        continue;
      }
      if (overlap >= n3->index - left)
      {
        continue;
      }
      if (n1->trace_back == -1)
      {
        continue;
      }
      if (overlap >= n3->stop_val - nodes[n1->trace_back].index - 2)
      {
        continue;
      }
      if ((stage == 1 && n3->cscore + n3->sscore +
          intergenic_mod(n3, n2, start_weight) > max_val) || (stage == 0 &&
          frame_bias[0]*n3->gc_score[0] + frame_bias[1]*n3->gc_score[1] +
          frame_bias[2]*n3->gc_score[2] > max_val))
      {
        max_frame = i;
        max_val = n3->cscore + n3->sscore +
                  intergenic_mod(n3, n2, start_weight);
      }
    }
    if (max_frame != -1)
    {
      n3 = &(nodes[n2->star_ptr[max_frame]]);
      if (stage == 0)
      {
        score = frame_bias[0]*n3->gc_score[0] + frame_bias[1]*n3->gc_score[1] +
                frame_bias[2]*n3->gc_score[2];
      }
      else if (stage == 1)
      {
        score = n3->cscore + n3->sscore + intergenic_mod(n3, n2, start_weight);
      }
    }
    else if (stage == 1)
    {
      score = intergenic_mod(n1, n2, start_weight);
    }
  }

  /* 5'rev->3'rev */
  else if (n1->strand == -1 && n1->type == START && n2->strand == -1 &&
           n2->type == STOP)
  {
    right -= 2;
    if (left >= right)
    {
      return;
    }
    if (stage == 1)
    {
      score = intergenic_mod(n1, n2, start_weight);
    }
  }

  /* 5'rev->5'fwd */
  else if (n1->strand == -1 && n1->type == START && n2->strand == 1 &&
           n2->type == START)
  {
    if (left >= right)
    {
      return;
    }
    if (stage == 1)
    {
      score = intergenic_mod(n1, n2, start_weight);
    }
  }

  /********************/
  /* Possible Operons */
  /********************/

  /* 3'fwd->3'fwd, check for a start just to left of first 3' */
  else if (n1->strand == 1 && n2->strand == 1 && n1->type == STOP &&
           n2->type == STOP)
  {
    if (n2->stop_val >= n1->index)
    {
      return;
    }
    if (n1->star_ptr[n2->index%3] == -1)
    {
      return;
    }
    n3 = &(nodes[n1->star_ptr[n2->index%3]]);
    left = n3->index;
    right += 2;
    if (stage == 0)
    {
      score = frame_bias[0]*n3->gc_score[0] + frame_bias[1]*n3->gc_score[1] +
              frame_bias[2]*n3->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n3->cscore + n3->sscore + intergenic_mod(n1, n3, start_weight);
    }
  }

  /* 3'rev->3'rev, check for a start just to right of second 3' */
  else if (n1->strand == -1 && n1->type == STOP && n2->strand == -1 &&
           n2->type == STOP)
  {
    if (n1->stop_val <= n2->index)
    {
      return;
    }
    if (n2->star_ptr[n1->index%3] == -1)
    {
      return;
    }
    n3 = &(nodes[n2->star_ptr[n1->index%3]]);
    left -= 2;
    right = n3->index;
    if (stage == 0)
    {
      score = frame_bias[0]*n3->gc_score[0] + frame_bias[1]*n3->gc_score[1] +
              frame_bias[2]*n3->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n3->cscore + n3->sscore + intergenic_mod(n3, n2, start_weight);
    }
  }

  /***************************************/
  /* Overlapping Opposite Strand 3' Ends */
  /***************************************/

  /* 3'for->5'rev */
  else if (n1->strand == 1 && n1->type == STOP && n2->strand == -1 &&
           n2->type == START)
  {
    if (n2->stop_val-2 >= n1->index+2)
    {
      return;
    }
    overlap = (n1->index+2) - (n2->stop_val-2) + 1;
    if (overlap >= MAX_OPP_OVLP)
    {
      return;
    }
    if ((n1->index+2 - n2->stop_val-2 + 1) >= (n2->index -n1->index+3 + 1))
    {
      return;
    }
    if (n1->trace_back == -1)
    {
      bound = 0;
    }
    else
    {
      bound = nodes[n1->trace_back].index;
    }
    if ((n1->index+2 - n2->stop_val-2 + 1) >= (n2->stop_val-3 - bound + 1))
    {
      return;
    }
    left = n2->stop_val-2;
    if (stage == 0)
    {
      score = frame_bias[0]*n2->gc_score[0] + frame_bias[1]*n2->gc_score[1] +
              frame_bias[2]*n2->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n2->cscore + n2->sscore - 0.15*start_weight;
    }
  }

  /* Frame bias score is multiplied by length */
  if (stage == 0)
  {
    score = ((double)(right-left+1-(overlap*2))) * score;
  }

  /* New best score for this node */
  if (n1->score + score >= n2->score)
  {
    n2->score = n1->score + score;
    n2->trace_back = node_a;
    n2->overlap_frame = max_frame;
  }

  return;
}

/******************************************************************************
  Dynamic programming is done but we now have to untangle the overlaps and
  connect nodes in the proper order (i.e. start to stop to start to stop),
  which sometimes involves going backwards due to overlapping genes.
******************************************************************************/
void untangle_overlaps(struct _node *nodes, int last_node)
{
  int i = 0;
  int path = last_node;
  int tmp = 0;
  int next_node = 0;

  if (path == -1)
  {
    return;
  }

  /* First Pass: untangle the triple overlaps, cases where two
     consecutive genes overlap with a third gene inside the
     overlap */
  while (nodes[path].trace_back != -1)
  {
    next_node = nodes[path].trace_back;
    if (nodes[path].strand == -1 && nodes[path].type == STOP &&
        nodes[next_node].strand == 1 && nodes[next_node].type == STOP &&
        nodes[path].overlap_frame != -1 &&
        nodes[path].index > nodes[next_node].index)
    {
      tmp = nodes[path].star_ptr[nodes[path].overlap_frame];
      for (i = tmp; nodes[i].index != nodes[tmp].stop_val; i--)
      {
      }
      nodes[path].trace_back = tmp;
      nodes[tmp].trace_back = i;
      nodes[i].overlap_frame = -1;
      nodes[i].trace_back = next_node;
    }
    path = nodes[path].trace_back;
  }

  /* Second Pass: Untangle the simple overlaps */
  path = last_node;
  while (nodes[path].trace_back != -1)
  {
    next_node = nodes[path].trace_back;
    if (nodes[path].strand == -1 && nodes[path].type == START &&
        nodes[next_node].strand == 1 && nodes[next_node].type == STOP)
    {
      for (i = path; nodes[i].index != nodes[path].stop_val; i--)
      {
      }
      nodes[path].trace_back = i;
      nodes[i].trace_back = next_node;
    }
    if (nodes[path].strand == 1 && nodes[path].type == STOP &&
        nodes[next_node].strand == 1 && nodes[next_node].type == STOP)
    {
      nodes[path].trace_back =
        nodes[next_node].star_ptr[(nodes[path].index)%3];
      nodes[nodes[path].trace_back].trace_back = next_node;
    }
    if (nodes[path].strand == -1 && nodes[path].type == STOP &&
        nodes[next_node].strand == -1 && nodes[next_node].type == STOP)
    {
      nodes[path].trace_back =
        nodes[path].star_ptr[(nodes[next_node].index)%3];
      nodes[nodes[path].trace_back].trace_back = next_node;
    }
    path = nodes[path].trace_back;
  }
}

/* Record trace forward pointers in the model */
void record_trace_forward_pointers(struct _node *nodes, int last_node)
{
  int path = last_node;

  if (path == -1)
  {
    return;
  }
  while (nodes[path].trace_back != -1)
  {
    nodes[nodes[path].trace_back].trace_forward = path;
    path = nodes[path].trace_back;
  }
}

/* Record trace forward pointers in the model */
void mark_valid_gene_nodes(struct _node *nodes, int last_node)
{
  int path = last_node;

  if (path == -1)
  {
    return;
  }
  nodes[path].status = 1;
  while (nodes[path].trace_back != -1)
  {
    path = nodes[path].trace_back;
    nodes[path].status = 1;
  }
}

/* Given the last node in the dynamic programming, find the */
/* first node using trace back. */
int find_first_node_from_last_node(struct _node *nodes, int last_node)
{
  int path = last_node;

  if (path == -1)
  {
    return -1;
  }
  while (nodes[path].trace_back != -1)
  {
    path = nodes[path].trace_back;
  }
  return path;
}
