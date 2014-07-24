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
  int i, j, initial_node, max_index = -1, path, nxt, tmp;
  double max_score = -1.0;

  if (num_nodes == 0)
  {
    return -1;
  }
  for (i = 0; i < num_nodes; i++)
  {
    nodes[i].score = 0;
    nodes[i].trace_back = -1;
    nodes[i].trace_forward = -1;
  }
  for (i = 0; i < num_nodes; i++)
  {

    /* Dynamic programming is n log n, so we set up a distance constraint */
    /* so that we will only try to connect two nodes within MAX_NODE_DIST */
    /* of each other in the sorted list of nodes.  We make an exception   */
    /* for same-gene connections that represent giant ORFs. */
    if (i < MAX_NODE_DIST)
    {
      initial_node = 0;
    }
    else
    {
      initial_node = i-MAX_NODE_DIST;
    }
    if (nodes[i].strand == -1 && is_start_node(&nodes[i]) == 1 &&
        nodes[initial_node].index >= nodes[i].stop_val)
    {
      while (initial_node >= 0 &&
             nodes[initial_node].index != nodes[i].stop_val)
      {
        initial_node--;
      }
    }
    if (nodes[i].strand == 1 && is_stop_node(&nodes[i]) == 1 &&
        nodes[initial_node].index >= nodes[i].stop_val)
    {
      while (initial_node >= 0 &&
             nodes[initial_node].index != nodes[i].stop_val)
      {
        initial_node--;
      }
    }
    if (initial_node < MAX_NODE_DIST)
    {
      initial_node = 0;
    }
    else
    {
      initial_node = initial_node-MAX_NODE_DIST;
    }
    /* Score all connections from nodes between initial_node and the current
       node-1 with the current_node */
    for (j = initial_node; j < i; j++)
    {
      score_connection(nodes, j, i, frame_bias, start_weight, stage);
    }
  }
  /* Locate the highest scoring node in the traceback and record its index */
  for (i = num_nodes-1; i >= 0; i--)
  {
    if (nodes[i].strand == 1 && is_start_node(&nodes[i]) == 1)
    {
      continue;
    }
    if (nodes[i].strand == -1 && is_stop_node(&nodes[i]) == 1)
    {
      continue;
    }
    if (nodes[i].score > max_score)
    {
      max_score = nodes[i].score;
      max_index = i;
    }
  }

  /* Dynamic programming is done but we now have to untangle   */
  /* the overlaps and connect nodes in the proper order (i.e.  */
  /* start to stop to start to stop), which sometimes involves */
  /* going backwards due to overlapping genes. */

  /* First Pass: untangle the triple overlaps, cases where two
     consecutive genes overlap with a third gene inside the
     overlap */
  path = max_index;
  while (nodes[path].trace_back != -1)
  {
    nxt = nodes[path].trace_back;
    if (nodes[path].strand == -1 && is_stop_node(&nodes[path]) == 1 &&
        nodes[nxt].strand == 1 && is_stop_node(&nodes[nxt]) == 1 &&
        nodes[path].overlap_frame != -1 &&
        nodes[path].index > nodes[nxt].index)
    {
      tmp = nodes[path].star_ptr[nodes[path].overlap_frame];
      for (i = tmp; nodes[i].index != nodes[tmp].stop_val; i--)
      {
      }
      nodes[path].trace_back = tmp;
      nodes[tmp].trace_back = i;
      nodes[i].overlap_frame = -1;
      nodes[i].trace_back = nxt;
    }
    path = nodes[path].trace_back;
  }

  /* Second Pass: Untangle the simple overlaps */
  path = max_index;
  while (nodes[path].trace_back != -1)
  {
    nxt = nodes[path].trace_back;
    if (nodes[path].strand == -1 && is_start_node(&nodes[path]) == 1 &&
        nodes[nxt].strand == 1 && is_stop_node(&nodes[nxt]) == 1)
    {
      for (i = path; nodes[i].index != nodes[path].stop_val; i--)
      {
      }
      nodes[path].trace_back = i;
      nodes[i].trace_back = nxt;
    }
    if (nodes[path].strand == 1 && is_stop_node(&nodes[path]) == 1 &&
        nodes[nxt].strand == 1 && is_stop_node(&nodes[nxt]) == 1)
    {
      nodes[path].trace_back = nodes[nxt].star_ptr[(nodes[path].index)%3];
      nodes[nodes[path].trace_back].trace_back = nxt;
    }
    if (nodes[path].strand == -1 && is_stop_node(&nodes[path]) == 1 &&
        nodes[nxt].strand == -1 && is_stop_node(&nodes[nxt]) == 1)
    {
      nodes[path].trace_back = nodes[path].star_ptr[(nodes[nxt].index)%3];
      nodes[nodes[path].trace_back].trace_back = nxt;
    }
    path = nodes[path].trace_back;
  }

  /* Dynamic programming did a trace_backack, but now we want to record
     trace-forward pointers, so we can get the genes in forward order. */
  path = max_index;
  while (nodes[path].trace_back != -1)
  {
    nodes[nodes[path].trace_back].trace_forward = path;
    path = nodes[path].trace_back;
  }

  /* Return the starting node of the dynamic programming or -1 if
     we found no genes in this sequence. */
  if (nodes[max_index].trace_back == -1)
  {
    return -1;
  }
  else
  {
    return max_index;
  }
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
  struct _node *n1 = &(nodes[node_a]);
  struct _node *n2 = &(nodes[node_b]);
  struct _node *n3 = NULL;
  int i, left = n1->index, right = n2->index, bound;
  int overlap = 0, max_frame = -1;
  double score = 0.0, score_mod = 0.0, max_val;

  /***********************/
  /* Invalid Connections */
  /***********************/

  /* 5'fwd->5'fwd, 5'rev->5'rev */
  if (is_start_node(n1) == 1 && is_start_node(n2) == 1 &&
      n1->strand == n2->strand)
  {
    return;
  }

  /* 5'fwd->5'rev, 5'fwd->3'rev */
  else if (n1->strand == 1 && is_start_node(n1) == 1 && n2->strand == -1)
  {
    return;
  }

  /* 3'rev->5'fwd, 3'rev->3'fwd) */
  else if (n1->strand == -1 && is_stop_node(n1) == 1 && n2->strand == 1)
  {
    return;
  }

  /* 5'rev->3'fwd */
  else if (n1->strand == -1 && is_start_node(n1) == 1 && n2->strand == 1 &&
           is_stop_node(n2) == 1)
  {
    return;
  }

  /******************/
  /* Edge Artifacts */
  /******************/
  if (n1->trace_back == -1 && n1->strand == 1 && is_stop_node(n1) == 1)
  {
    return;
  }
  if (n1->trace_back == -1 && n1->strand == -1 && is_start_node(n1) == 1)
  {
    return;
  }

  /*********/
  /* Genes */
  /*********/

  /* 5'fwd->3'fwd */
  else if (n1->strand == n2->strand && n1->strand == 1 &&
           is_start_node(n1) == 1 && is_stop_node(n2) == 1)
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
      score_mod = frame_bias[0]*n1->gc_score[0] +
                  frame_bias[1]*n1->gc_score[1] +
                  frame_bias[2]*n1->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n1->cscore + n1->sscore;
    }
  }

  /* 3'rev->5'rev */
  else if (n1->strand == n2->strand && n1->strand == -1 &&
           is_stop_node(n1) == 1 && is_start_node(n2) == 1)
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
      score_mod = frame_bias[0]*n2->gc_score[0] +
                  frame_bias[1]*n2->gc_score[1] +
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
  else if (n1->strand == 1 && is_stop_node(n1) == 1 && n2->strand == 1 &&
           is_start_node(n2) == 1)
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
  else if (n1->strand == 1 && is_stop_node(n1) == 1 && n2->strand == -1 &&
           is_stop_node(n2) == 1)
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
        score_mod = frame_bias[0]*n3->gc_score[0] +
                    frame_bias[1]*n3->gc_score[1] +
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
  else if (n1->strand == -1 && is_start_node(n1) == 1 && n2->strand == -1 &&
           is_stop_node(n2) == 1)
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
  else if (n1->strand == -1 && is_start_node(n1) == 1 && n2->strand == 1 &&
           is_start_node(n2) == 1)
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
  else if (n1->strand == 1 && n2->strand == 1 && is_stop_node(n1) == 1 &&
           is_stop_node(n2) == 1)
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
      score_mod = frame_bias[0]*n3->gc_score[0] +
                  frame_bias[1]*n3->gc_score[1] +
                  frame_bias[2]*n3->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n3->cscore + n3->sscore + intergenic_mod(n1, n3, start_weight);
    }
  }

  /* 3'rev->3'rev, check for a start just to right of second 3' */
  else if (n1->strand == -1 && is_stop_node(n1) == 1 && n2->strand == -1 &&
           is_stop_node(n2) == 1)
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
      score_mod = frame_bias[0]*n3->gc_score[0] +
                  frame_bias[1]*n3->gc_score[1] +
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
  else if (n1->strand == 1 && is_stop_node(n1) == 1 && n2->strand == -1 &&
           is_start_node(n2) == 1)
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
      score_mod = frame_bias[0]*n2->gc_score[0] +
                  frame_bias[1]*n2->gc_score[1] +
                  frame_bias[2]*n2->gc_score[2];
    }
    else if (stage == 1)
    {
      score = n2->cscore + n2->sscore - 0.15*start_weight;
    }
  }

  if (stage == 0)
  {
    score = ((double)(right-left+1-(overlap*2)))*score_mod;
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
  Sometimes bad genes creep into the model due to the node distance constraint
  in the dynamic programming routine.  This routine just does a sweep through
  the genes and eliminates ones with negative scores.  The "elim" flag is
  set to 1 for nodes that should be ignored when building the list of genes.
******************************************************************************/
void eliminate_bad_genes(struct _node *nodes, int initial_node,
                         double start_weight)
{
  int path;

  if (initial_node == -1)
  {
    return;
  }
  path = initial_node;
  while (nodes[path].trace_back != -1)
  {
    path = nodes[path].trace_back;
  }
  while (nodes[path].trace_forward != -1)
  {
    if (nodes[path].strand == 1 && is_stop_node(&nodes[path]) == 1)
    {
      nodes[nodes[path].trace_forward].sscore +=
        intergenic_mod(&nodes[path], &nodes[nodes[path].trace_forward],
                       start_weight);
    }
    if (nodes[path].strand == -1 && is_start_node(&nodes[path]) == 1)
    {
      nodes[path].sscore +=
        intergenic_mod(&nodes[path], &nodes[nodes[path].trace_forward],
                       start_weight);
    }
    path = nodes[path].trace_forward;
  }

  path = initial_node;
  while (nodes[path].trace_back != -1)
  {
    path = nodes[path].trace_back;
  }
  while (nodes[path].trace_forward != -1)
  {
    if (nodes[path].strand == 1 && is_start_node(&nodes[path]) == 1 &&
        nodes[path].cscore + nodes[path].sscore < 0)
    {
      nodes[path].eliminate = 1;
      nodes[nodes[path].trace_forward].eliminate = 1;
    }
    if (nodes[path].strand == -1 && is_stop_node(&nodes[path]) == 1 &&
        nodes[nodes[path].trace_forward].cscore +
        nodes[nodes[path].trace_forward].sscore < 0)
    {
      nodes[path].eliminate = 1;
      nodes[nodes[path].trace_forward].eliminate = 1;
    }
    path = nodes[path].trace_forward;
  }
}
