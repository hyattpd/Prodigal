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

#include "node.h"

/*******************************************************************************
  Adds nodes to the node list.  Genes must be >=90bp in length, unless they
  run off the edge, in which case they only have to be 60bp.
*******************************************************************************/

int add_nodes(unsigned char *seq, unsigned char *rseq, unsigned char *useq, 
              int slen, struct _node *nodes, int closed, int cross_gaps, 
              int trans_table) {
  int i, nn = 0, last[3], saw_start[3], min_dist[3], edge[3];
  int slmod = 0;

  /* Forward strand nodes */
  slmod = slen%3;
  for(i = 0; i < 3; i++) {
    last[(i+slmod)%3] = slen+i; 
    while(last[(i+slmod)%3]+2 > slen-1) last[(i+slmod)%3]-=3;
  }
  for(i = 0; i < 3; i++) {
    saw_start[i%3] = 0;
    if(is_stop(seq, last[i%3], trans_table) == 1) {
      min_dist[i%3] = MIN_GENE;
      edge[i%3] = 0;
    }
    else {
      min_dist[i%3] = MIN_EDGE_GENE;
      edge[i%3] = 1;
    }
  }
  for(i = slen-3; i >= 0; i--) {
    if(is_stop(seq, i, trans_table) == 1 || (i <= slen-12 && cross_gaps == 0 &&
       gap_to_right(useq, i) == 1)) {
      if(saw_start[i%3] == 1 && (closed == 0 || edge[i%3] == 0)) {
        nodes[nn].edge = edge[i%3];
        nodes[nn].ndx = last[i%3]; 
        nodes[nn].type = STOP;
        nodes[nn].strand = 1; 
        nodes[nn++].stop_val = i;
      }
      if(is_stop(seq, i, trans_table) == 1) {
        min_dist[i%3] = MIN_GENE;
        edge[i%3] = 0;
      }
      else {
        min_dist[i%3] = MIN_EDGE_GENE;
        edge[i%3] = 1;
      }
      last[i%3] = i; 
      saw_start[i%3] = 0;
      continue;
    }
    if(edge[i%3] == 1 && closed == 1) continue;
     
    if(is_start(seq, i, trans_table) == 1 && 
       ((last[i%3]-i+3) >= min_dist[i%3])) {
      nodes[nn].ndx = i; 
      if(is_atg(seq, i) == 1) nodes[nn].type = ATG; 
      else if(is_gtg(seq, i) == 1) nodes[nn].type = GTG; 
      else if(is_ttg(seq, i) == 1) nodes[nn].type = TTG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = last[i%3]; 
      nodes[nn++].strand = 1;
    }
    else if(closed == 0 && ((last[i%3]-i) > MIN_EDGE_GENE) && (i <= 2 ||
            (cross_gaps == 0 && i >= 9 && codon_has_n(useq, i) == 0 &&
            gap_to_left(useq, i) == 1))) {
      nodes[nn].ndx = i; 
      nodes[nn].type = ATG; 
      saw_start[i%3] = 1;
      nodes[nn].edge = 1; 
      nodes[nn].stop_val = last[i%3];
      nodes[nn++].strand = 1;
    }
  }
  for(i = 0; i < 3; i++) {
    if(saw_start[i%3] == 1) {
      nodes[nn].edge = edge[i%3];
      nodes[nn].ndx = last[i%3]; 
      nodes[nn].type = STOP;
      nodes[nn].strand = 1; 
      nodes[nn++].stop_val = i-6;
    }
  }

  /* Reverse strand nodes */
  for(i = 0; i < 3; i++) {
    last[(i+slmod)%3] = slen+i; 
    while(last[(i+slmod)%3]+2 > slen-1) last[(i+slmod)%3]-=3;
  }
  for(i = 0; i < 3; i++) {
    saw_start[i%3] = 0;
    if(is_stop(rseq, last[i%3], trans_table) == 1) {
      min_dist[i%3] = MIN_GENE;
      edge[i%3] = 1;
    } 
    else {
      min_dist[i%3] = MIN_EDGE_GENE;
      edge[i%3] = 1;
    }
  }
  for(i = slen-3; i >= 0; i--) {
    if(is_stop(rseq, i, trans_table)==1 || (i <= slen-12 && cross_gaps == 0 &&
       gap_to_left(useq, slen-3-i) == 1)) {
      if(saw_start[i%3] == 1 && (closed == 0 || edge[i%3] == 0)) {
        nodes[nn].edge = edge[i%3];
        nodes[nn].ndx = slen-last[i%3]-1; 
        nodes[nn].type = STOP;
        nodes[nn].strand = -1; 
        nodes[nn++].stop_val = slen-i-1;
      }
      if(is_stop(rseq, i, trans_table) == 1) {
        min_dist[i%3] = MIN_GENE;
        edge[i%3] = 0;
      }
      else {
        min_dist[i%3] = MIN_EDGE_GENE;
        edge[i%3] = 1;
      }
      last[i%3] = i; 
      saw_start[i%3] = 0;
      continue;
    }
    if(edge[i%3] == 1 && closed == 1) continue;

    if(is_start(rseq, i, trans_table) == 1 && 
       ((last[i%3]-i+3) >= min_dist[i%3])) {
      nodes[nn].ndx = slen-i-1; 
      if(is_atg(rseq, i) == 1) nodes[nn].type = ATG; 
      else if(is_gtg(rseq, i) == 1) nodes[nn].type = GTG; 
      else if(is_ttg(rseq, i) == 1) nodes[nn].type = TTG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = slen-last[i%3]-1; 
      nodes[nn++].strand = -1;
    }
    else if(closed == 0 && ((last[i%3]-i) > MIN_EDGE_GENE) && (i <= 2 ||
            (cross_gaps == 0 && i >= 9 && codon_has_n(useq, slen-3-i) == 0 &&
            gap_to_right(useq, slen-3-i) == 1))) {
      nodes[nn].ndx = slen-i-1; 
      nodes[nn].type = ATG; 
      saw_start[i%3] = 1;
      nodes[nn].edge = 1; 
      nodes[nn].stop_val = slen-last[i%3]-1;
      nodes[nn++].strand = -1;
    }
  }
  for(i = 0; i < 3; i++) {
    if(saw_start[i%3] == 1) {
      nodes[nn].edge = edge[i%3];
      nodes[nn].ndx = slen - last[i%3] - 1; 
      nodes[nn].type = STOP;
      nodes[nn].strand = -1; 
      nodes[nn++].stop_val = slen-i+5;
    }
  }
  return nn;
}

/* Memset nodes to 0 and return 0 */
int zero_nodes(struct _node *nod, int nn) {
  memset(nod, 0, nn*sizeof(struct _node));
  return 0;
}

/* Simple routine to zero out the node scores */

void reset_node_scores(struct _node *nod, int nn) {
  int i, j;
  for(i = 0; i < nn; i++) {
    for(j = 0; j < 3; j++) {
      nod[i].star_ptr[j] = 0;
      nod[i].gc_score[j] = 0.0;
    }
    for(j = 0; j < 2; j++) nod[i].rbs[j] = 0;
    nod[i].score = 0.0;
    nod[i].cscore = 0.0; 
    nod[i].sscore = 0.0; 
    nod[i].rscore = 0.0; 
    nod[i].tscore = 0.0; 
    nod[i].uscore = 0.0; 
    nod[i].traceb = -1;
    nod[i].tracef = -1;
    nod[i].ov_mark = -1;
    nod[i].elim = 0;
    nod[i].gc_bias = 0;
    memset(&nod[i].mot, 0, sizeof(struct _motif));
  }
}

/*******************************************************************************
  Since dynamic programming can't go 'backwards', we have to record
  information about overlapping genes in order to build the models.  So, for
  example, in cases like 5'->3', 5'-3' overlapping on the same strand, we
  record information about the 2nd 5' end under the first 3' end's
  information.  For every stop, we calculate and store all the best starts
  that could be used in genes that overlap that 3' end.
*******************************************************************************/

void record_overlapping_starts(struct _node *nod, int nn, double st_wt, int
                               flag) {
  int i, j;
  double max_sc;

  for(i = 0; i < nn; i++) {
    for(j = 0; j < 3; j++) nod[i].star_ptr[j] = -1;
    if(nod[i].type != STOP || nod[i].edge == 1) continue;
    if(nod[i].strand == 1) {
      max_sc = -100;
      for(j = i+3; j >= 0; j--) {
        if(j >= nn || nod[j].ndx > nod[i].ndx+2) continue;
        if(nod[j].ndx + MAX_SAM_OVLP < nod[i].ndx) break;
        if(nod[j].strand == 1 && nod[j].type != STOP) {
          if(nod[j].stop_val <= nod[i].ndx) continue; 
          if(flag == 0 && nod[i].star_ptr[(nod[j].ndx)%3] == -1)
            nod[i].star_ptr[(nod[j].ndx)%3] = j;
          else if(flag == 1 && (nod[j].cscore + nod[j].sscore +
                  intergenic_mod(&nod[i], &nod[j], st_wt) > max_sc)) {
            nod[i].star_ptr[(nod[j].ndx)%3] = j;
            max_sc = nod[j].cscore + nod[j].sscore +
                     intergenic_mod(&nod[i], &nod[j], st_wt);
          }
        }
      }
    }
    else {
      max_sc = -100;
      for(j = i-3; j < nn; j++) {
        if(j < 0 || nod[j].ndx < nod[i].ndx-2) continue;
        if(nod[j].ndx - MAX_SAM_OVLP > nod[i].ndx) break;
        if(nod[j].strand == -1 && nod[j].type != STOP) {
          if(nod[j].stop_val >= nod[i].ndx) continue; 
          if(flag == 0 && nod[i].star_ptr[(nod[j].ndx)%3] == -1)
            nod[i].star_ptr[(nod[j].ndx)%3] = j;
          else if(flag == 1 && (nod[j].cscore + nod[j].sscore +
                  intergenic_mod(&nod[j], &nod[i], st_wt) > max_sc)) {
            nod[i].star_ptr[(nod[j].ndx)%3] = j;
            max_sc = nod[j].cscore + nod[j].sscore +
                     intergenic_mod(&nod[j], &nod[i], st_wt);
          }
        }
      }
    }
  }
}

/*******************************************************************************
  This routine goes through all the ORFs and counts the relative frequency of
  the most common frame for G+C content.  In high GC genomes, this tends to be
  the third position.  In low GC genomes, this tends to be the first position.
  Genes will be selected as a training set based on the nature of this bias
  for this particular organism.
*******************************************************************************/

void frame_score(int *gc, struct _node *nod, int nn) {
  int i, j, ctr[3][3], last[3], frmod, fr, mfr;

  if(nn == 0) return;
  for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) ctr[i][j] = 0;
  for(i = nn-1; i >= 0; i--) {
    fr = (nod[i].ndx)%3; frmod = 3 - fr;
    if(nod[i].strand == 1 && nod[i].type == STOP) {
      for(j = 0; j < 3; j++) ctr[fr][j] = 0;
      last[fr] = nod[i].ndx;
      ctr[fr][(gc[nod[i].ndx] + frmod)%3] = 1;
    }
    else if(nod[i].strand == 1) {
      for(j = last[fr]-3; j >= nod[i].ndx; j-=3) ctr[fr][(gc[j] + frmod)%3] ++;
      mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]);
      nod[i].gc_bias = mfr;
      for(j = 0; j < 3; j++) {
        nod[i].gc_score[j] = (3.0*ctr[fr][j]);
        nod[i].gc_score[j] /= (1.0*(nod[i].stop_val - nod[i].ndx + 3));
      }
      last[fr] = nod[i].ndx;
    }
  }
  for(i = 0; i < nn; i++) {
    fr = (nod[i].ndx)%3; frmod = fr;
    if(nod[i].strand == -1 && nod[i].type == STOP) {
      for(j = 0; j < 3; j++) ctr[fr][j] = 0;
      last[fr] = nod[i].ndx;
      ctr[fr][((3-gc[nod[i].ndx]) + frmod)%3] = 1;
    }
    else if(nod[i].strand == -1) {
      for(j = last[fr]+3; j <= nod[i].ndx; j+=3)
        ctr[fr][((3-gc[j]) + frmod)%3]++;
      mfr = max_fr(ctr[fr][0], ctr[fr][1], ctr[fr][2]);
      nod[i].gc_bias = mfr;
      for(j = 0; j < 3; j++) {
        nod[i].gc_score[j] = (3.0*ctr[fr][j]);
        nod[i].gc_score[j] /= (1.0*(nod[i].ndx - nod[i].stop_val + 3));
      }
      last[fr] = nod[i].ndx;
    }
  }
}

/*******************************************************************************
  Scoring function for all the start nodes.  This score has two factors:  (1)
  Coding, which is a composite of coding score and length, and (2) Start
  score, which is a composite of RBS score and ATG/TTG/GTG.
*******************************************************************************/

void score_nodes(unsigned char *seq, unsigned char *rseq, int slen,
                 struct _node *nod, int nn, struct _training *tinf,
                 int closed, int mode) {
  int i, j;
  double negf, posf, rbs1, rbs2, sd_score, edge_gene, min_anon_len;

  /* Step 1: Calculate raw coding potential for every start-stop pair. */
  calc_orf_gc(seq, rseq, slen, nod, nn);
  raw_coding_score(seq, rseq, slen, nod, nn, tinf->trans_table, tinf->gc, tinf->gene_dc);

  /* Step 2: Calculate raw RBS Scores for every start node. */
  if(tinf->uses_sd == 1) rbs_score(seq, rseq, slen, nod, nn, tinf->rbs_wt);
  else {
    for(i = 0; i < nn; i++) {
      if(nod[i].type == STOP || nod[i].edge == 1) continue;
      find_best_upstream_motif(tinf, seq, rseq, slen, &nod[i], 2);
    }
  }

  /* Step 3: Score the start nodes */
  for(i = 0; i < nn; i++) {
    if(nod[i].type == STOP) continue;

    /* Does this gene run off the edge? */
    edge_gene = 0;
    if(nod[i].edge == 1) edge_gene++;
    if((nod[i].strand == 1 && is_stop(seq, nod[i].stop_val,
       tinf->trans_table) == 0) || (nod[i].strand == -1 && is_stop(rseq, slen-1-
       nod[i].stop_val, tinf->trans_table) == 0)) edge_gene++;

    /* Edge Nodes : stops with no starts, give a small bonus */
    if(nod[i].edge == 1) {
      nod[i].tscore = EDGE_BONUS*tinf->st_wt/edge_gene;
      nod[i].uscore = 0.0;
      nod[i].rscore = 0.0;
    }

    else {

      /* Type Score */
      nod[i].tscore = tinf->type_wt[nod[i].type] * tinf->st_wt;

      /* RBS Motif Score */
      rbs1 = tinf->rbs_wt[nod[i].rbs[0]];
      rbs2 = tinf->rbs_wt[nod[i].rbs[1]];
      sd_score = dmax(rbs1, rbs2) * tinf->st_wt;
      if(tinf->uses_sd == 1) nod[i].rscore = sd_score;
      else {
        nod[i].rscore = tinf->st_wt*nod[i].mot.score;
        if(nod[i].rscore < sd_score && tinf->no_mot > -0.5)
          nod[i].rscore = sd_score;
      }

      /* Upstream Score */
      if(nod[i].strand == 1) 
        score_upstream_composition(seq, slen, &nod[i], tinf);
      else score_upstream_composition(rseq, slen, &nod[i], tinf);

      /****************************************************************
      ** Penalize upstream score if choosing this start would stop   **
      ** the gene from running off the edge.                         **
      ****************************************************************/
      if(closed == 0 && nod[i].ndx <= 2 && nod[i].strand == 1) 
        nod[i].uscore += EDGE_UPS*tinf->st_wt; 
      else if(closed == 0 && nod[i].ndx >= slen-3 && nod[i].strand == -1)
        nod[i].uscore += EDGE_UPS*tinf->st_wt; 
      else if(i < 500 && nod[i].strand == 1) {
        for(j = i-1; j >= 0; j--)
          if(nod[j].edge == 1 && nod[i].stop_val == nod[j].stop_val) {
            nod[i].uscore += EDGE_UPS*tinf->st_wt; 
            break;
          }
      }
      else if(i >= nn-500 && nod[i].strand == -1) {
        for(j = i+1; j < nn; j++)
          if(nod[j].edge == 1 && nod[i].stop_val == nod[j].stop_val) {
            nod[i].uscore += EDGE_UPS*tinf->st_wt; 
            break;
          }
      }

    }

    /* Convert starts at base 1 and slen to edge genes if closed = 0 */
    if(((nod[i].ndx <= 2 && nod[i].strand == 1) || (nod[i].ndx >= slen-3 &&
       nod[i].strand == -1)) && nod[i].edge == 0 && closed == 0) {
      edge_gene++;
      nod[i].edge = 1;
      nod[i].tscore = 0.0;
      nod[i].uscore = EDGE_BONUS*tinf->st_wt/edge_gene;
      nod[i].rscore = 0.0;
    }

    /* Penalize starts with no stop codon */
    if(nod[i].edge == 0 && edge_gene == 1) 
      nod[i].uscore -= 0.5*EDGE_BONUS*tinf->st_wt;

    /* Penalize non-edge genes < 250bp */
    if(edge_gene == 0 && abs(nod[i].ndx-nod[i].stop_val) < 250) {
      negf = 250.0/(float)abs(nod[i].ndx-nod[i].stop_val);
      posf = (float)abs(nod[i].ndx-nod[i].stop_val)/250.0;
      if(nod[i].rscore < 0) nod[i].rscore *= negf; 
      if(nod[i].uscore < 0) nod[i].uscore *= negf; 
      if(nod[i].tscore < 0) nod[i].tscore *= negf; 
      if(nod[i].rscore > 0) nod[i].rscore *= posf; 
      if(nod[i].uscore > 0) nod[i].uscore *= posf; 
      if(nod[i].tscore > 0) nod[i].tscore *= posf; 
    }

    /**************************************************************/
    /* Coding Penalization in Anonymous Fragments:    Internal    */
    /* genes must have a score of 5.0 and be >= 120bp.  High GC   */
    /* genes are also penalized.                                  */
    /**************************************************************/
    if(mode == MODE_ANON && slen < 3000 && edge_gene == 0 && 
       (nod[i].cscore < 5.0 || abs(nod[i].ndx-nod[i].stop_val < 120))) {
      nod[i].cscore -= META_PEN*dmax(0, (3000-slen)/2700.0);
    }
 
    /* Base Start Score */
    nod[i].sscore = nod[i].tscore + nod[i].rscore + nod[i].uscore;

    /**************************************************************/
    /* Penalize starts if coding is negative.  Larger penalty for */
    /* edge genes, since the start is offset by a smaller amount  */
    /* of coding than normal.                                     */
    /**************************************************************/
    if(nod[i].cscore < 0.0) {
      if(edge_gene > 0 && nod[i].edge == 0) {
        if(mode != MODE_ANON || slen > 1500) nod[i].sscore -= tinf->st_wt;
        else nod[i].sscore -= (10.31 - 0.004*slen);
      }
      else if(mode == MODE_ANON && slen < 3000 && nod[i].edge == 1) {
        min_anon_len = sqrt(slen)*5.0;
        if(abs(nod[i].ndx-nod[i].stop_val) >= min_anon_len) {
          if(nod[i].cscore >= 0) nod[i].cscore = -1.0;
          nod[i].sscore = 0.0; 
          nod[i].uscore = 0.0; 
        }
      }
      else nod[i].sscore -= 0.5;
    }
    else if(nod[i].cscore < 5.0 && mode == MODE_ANON && abs(nod[i].ndx-
            nod[i].stop_val < 120) && nod[i].sscore < 0.0)
      nod[i].sscore -= tinf->st_wt; 
  }
}

/* Calculate the GC Content for each start-stop pair */
void calc_orf_gc(unsigned char *seq, unsigned char *rseq, int slen, struct
                 _node *nod, int nn) {
  int i, j, last[3], fr;
  double gc[3], gsize = 0.0;

  /* Go through each start-stop pair and calculate the %GC of the gene */
  for(i = 0; i < 3; i++) gc[i] = 0.0;
  for(i = nn-1; i >= 0; i--) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == 1 && nod[i].type == STOP) {
      last[fr] = nod[i].ndx;
      gc[fr] = is_gc(seq, nod[i].ndx) + is_gc(seq, nod[i].ndx+1) +
               is_gc(seq, nod[i].ndx+2);
    }
    else if(nod[i].strand == 1) {
      for(j = last[fr]-3; j >= nod[i].ndx; j-=3)
        gc[fr] += is_gc(seq, j) + is_gc(seq, j+1) + is_gc(seq, j+2);
      gsize = (float)(abs(nod[i].stop_val-nod[i].ndx)+3.0);
      nod[i].gc_cont = gc[fr]/gsize;
      last[fr] = nod[i].ndx;
    }
  }
  for(i = 0; i < 3; i++) gc[i] = 0.0;
  for(i = 0; i < nn; i++) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == -1 && nod[i].type == STOP) {
      last[fr] = nod[i].ndx;
      gc[fr] = is_gc(seq, nod[i].ndx) + is_gc(seq, nod[i].ndx-1) +
               is_gc(seq, nod[i].ndx-2);
    }
    else if(nod[i].strand == -1) {
      for(j = last[fr]+3; j <= nod[i].ndx; j+=3)
        gc[fr] += is_gc(seq, j) + is_gc(seq, j+1) + is_gc(seq, j+2);
      gsize = (float)(abs(nod[i].stop_val-nod[i].ndx)+3.0);
      nod[i].gc_cont = gc[fr]/gsize;
      last[fr] = nod[i].ndx;
    }
  }
}

/*******************************************************************************
  Score each candidate's coding.  We also sharpen coding/noncoding thresholds
  to prevent choosing interior starts when there is strong coding continuing
  upstream.
*******************************************************************************/

void raw_coding_score(unsigned char *seq, unsigned char *rseq, int slen, struct
                      _node *nod, int nn, int trans_table, double gc, 
                      double *gene_dc) {
  int i, j, last[3], fr;
  double score[3], lfac, no_stop, gsize = 0.0;

  no_stop = 1.0-prob_stop(trans_table, gc);

  /* Initial Pass: Score coding potential (start->stop) */
  for(i = 0; i < 3; i++) score[i] = 0.0;
  for(i = nn-1; i >= 0; i--) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == 1 && nod[i].type == STOP) {
      last[fr] = nod[i].ndx;
      score[fr] = 0.0;
    }
    else if(nod[i].strand == 1) {
      for(j = last[fr]-3; j >= nod[i].ndx; j-=3)
        score[fr] += gene_dc[mer_ndx(6, seq, j)];
      nod[i].cscore = score[fr];
      last[fr] = nod[i].ndx;
    }
  }
  for(i = 0; i < 3; i++) score[i] = 0.0;
  for(i = 0; i < nn; i++) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == -1 && nod[i].type == STOP) {
      last[fr] = nod[i].ndx;
      score[fr] = 0.0;
    }
    else if(nod[i].strand == -1) {
      for(j = last[fr]+3; j <= nod[i].ndx; j+=3)
        score[fr] += gene_dc[mer_ndx(6, rseq, slen-j-1)];
      nod[i].cscore = score[fr];
      last[fr] = nod[i].ndx;
    }
  }

  /* Second Pass: Penalize start nodes with ascending coding to their left */
  for(i = 0; i < 3; i++) score[i] = -10000.0;
  for(i = 0; i < nn; i++) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == 1 && nod[i].type == STOP) score[fr] = -10000.0;
    else if(nod[i].strand == 1) {
      if(nod[i].cscore > score[fr]) score[fr] = nod[i].cscore;
      else nod[i].cscore -= (score[fr] - nod[i].cscore);
    }
  }
  for(i = 0; i < 3; i++) score[i] = -10000.0;
  for(i = nn-1; i >= 0; i--) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == -1 && nod[i].type == STOP) score[fr] = -10000.0;
    else if(nod[i].strand == -1) {
      if(nod[i].cscore > score[fr]) score[fr] = nod[i].cscore;
      else nod[i].cscore -= (score[fr] - nod[i].cscore);
    }
  }

  /* Third Pass: Add length-based factor to the score      */
  /* Penalize start nodes based on length to their left    */
  for(i = 0; i < nn; i++) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == 1 && nod[i].type == STOP) score[fr] = -10000.0;
    else if(nod[i].strand == 1) {
      gsize = ((float)(abs(nod[i].stop_val-nod[i].ndx)+3.0))/3.0;
      if(gsize > 1000.0) {
        lfac = log((1-pow(no_stop, 1000.0))/pow(no_stop, 1000.0));
        lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80));
        lfac *= (gsize - 80) / 920.0;
      }
      else {
        lfac = log((1-pow(no_stop, gsize))/pow(no_stop, gsize));
        lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80));
      }
      if(lfac > score[fr]) score[fr] = lfac;
      else lfac -= dmax(dmin(score[fr] - lfac, lfac), 0);
      if(lfac > 3.0 && nod[i].cscore < 0.5*lfac) nod[i].cscore = 0.5*lfac;
      nod[i].cscore += lfac;
    }
  }
  for(i = nn-1; i >= 0; i--) {
    fr = (nod[i].ndx)%3;
    if(nod[i].strand == -1 && nod[i].type == STOP) score[fr] = -10000.0;
    else if(nod[i].strand == -1) {
      gsize = ((float)(abs(nod[i].stop_val-nod[i].ndx)+3.0))/3.0;
      if(gsize > 1000.0) {
        lfac = log((1-pow(no_stop, 1000.0))/pow(no_stop, 1000.0));
        lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80));
        lfac *= (gsize - 80) / 920.0;
      }
      else {
        lfac = log((1-pow(no_stop, gsize))/pow(no_stop, gsize));
        lfac -= log((1-pow(no_stop, 80))/pow(no_stop, 80));
      }
      if(lfac > score[fr]) score[fr] = lfac;
      else lfac -= dmax(dmin(score[fr] - lfac, lfac), 0);
      if(lfac > 3.0 && nod[i].cscore < 0.5*lfac) nod[i].cscore = 0.5*lfac;
      nod[i].cscore += lfac;
    }
  }
}

/*******************************************************************************
  RBS Scoring Function: Calculate the RBS motif and then multiply it by the
  appropriate weight for that motif (determined in the start training
  function).
*******************************************************************************/
void rbs_score(unsigned char *seq, unsigned char *rseq, int slen, struct _node
               *nod, int nn, double *rbs_wt) {
  int i, j;
  int cur_sc[2];

  /* Scan all starts looking for RBS's */
  for(i = 0; i < nn; i++) {
    if(nod[i].type == STOP || nod[i].edge == 1) continue;
    nod[i].rbs[0] = 0;
    nod[i].rbs[1] = 0;
    if(nod[i].strand == 1) {
      for(j = nod[i].ndx - 20; j <= nod[i].ndx - 6; j++) {
        if(j < 0) continue;
        cur_sc[0] = shine_dalgarno_exact(seq, j, nod[i].ndx, rbs_wt);
        cur_sc[1] = shine_dalgarno_mm(seq, j, nod[i].ndx, rbs_wt);
        if(cur_sc[0] > nod[i].rbs[0]) nod[i].rbs[0] = cur_sc[0];
        if(cur_sc[1] > nod[i].rbs[1]) nod[i].rbs[1] = cur_sc[1];
      }
    }
    else if(nod[i].strand == -1) {
      for(j = slen - nod[i].ndx - 21; j <= slen - nod[i].ndx - 7; j++) {
        if(j > slen-1) continue;
        cur_sc[0] = shine_dalgarno_exact(rseq, j, slen-1-nod[i].ndx,
                                         rbs_wt);
        cur_sc[1] = shine_dalgarno_mm(rseq, j, slen-1-nod[i].ndx,
                                      rbs_wt);
        if(cur_sc[0] > nod[i].rbs[0]) nod[i].rbs[0] = cur_sc[0];
        if(cur_sc[1] > nod[i].rbs[1]) nod[i].rbs[1] = cur_sc[1];
      }
    }
  }
}

/*******************************************************************************
  For a given start, score the base composition of the upstream region at
  positions -1 and -2 and -15 to -44.  This will be used to supplement the
  SD (or other) motif finder with additional information.
*******************************************************************************/
void score_upstream_composition(unsigned char *seq, int slen, struct _node *nod,
                                struct _training *tinf) {
  int i, start, count = 0;
  if(nod->strand == 1) start = nod->ndx;
  else start = slen-1-nod->ndx;

  nod->uscore = 0.0;
  for(i = 1; i < 45; i++) {
    if(i > 2 && i < 15) continue;
    if(start-i < 0) continue;
    nod->uscore += 0.4*tinf->st_wt*
                   tinf->ups_comp[count][mer_ndx(1, seq, start-i)];
    count++;
  }
}

/*******************************************************************************
  Given the weights for various motifs/distances from the training file,
  return the highest scoring mer/spacer combination of 3-6bp motifs with a
  spacer ranging from 3bp to 15bp.  In the final stage of start training, only
  good scoring motifs are returned.
*******************************************************************************/
void find_best_upstream_motif(struct _training *tinf, unsigned char *seq,
                              unsigned char *rseq, int slen, struct _node *nod,
                              int stage) {
  int i, j, start, spacer, spacendx, index;
  int max_spacer = 0, max_spacendx = 0, max_len = 0, max_ndx = 0;
  double max_sc = -100.0, score = 0.0;
  unsigned char *wseq;

  if(nod->type == STOP || nod->edge == 1) return;
  if(nod->strand == 1) { wseq = seq; start = nod->ndx; }
  else { wseq = rseq; start = slen-1-nod->ndx; }

  for(i = 3; i >= 0; i--) {
    for(j = start-18-i; j <= start-6-i; j++) {
      if(j < 0) continue;
      spacer = start-j-i-3;
      if(j <= start-16-i) spacendx = 3;
      else if(j <= start-14-i) spacendx = 2;
      else if(j >= start-7-i) spacendx = 1;
      else spacendx = 0;
      index = mer_ndx(i+3, wseq, j);
      score = tinf->mot_wt[i][spacendx][index];
      if(score > max_sc) {
        max_sc = score;
        max_spacendx = spacendx;
        max_spacer = spacer;
        max_ndx = index;
        max_len = i+3;
      }
    }
  }

  if(stage == 2 && (max_sc == -4.0 || max_sc < tinf->no_mot + 0.69)) {
    nod->mot.ndx = 0;
    nod->mot.len = 0;
    nod->mot.spacendx = 0;
    nod->mot.spacer = 0;
    nod->mot.score = tinf->no_mot;
  }
  else {
    nod->mot.ndx = max_ndx;
    nod->mot.len = max_len;
    nod->mot.spacendx = max_spacendx;
    nod->mot.spacer = max_spacer;
    nod->mot.score = max_sc;
  }
}

/*******************************************************************************
  When connecting two genes, we add a bonus for the -1 and -4 base overlaps on
  the same strand, which often signify an operon and negate the need for an
  RBS for the second gene.  In addition, we add a slight bonus when genes are
  close and a slight penalty when switching strands or having a large
  intergenic space.
*******************************************************************************/
double intergenic_mod(struct _node *n1, struct _node *n2, double start_weight) {
  int dist;
  double rval = 0.0, ovlp = 0.0;
  if((n1->strand == 1 && n2->strand == 1 &&
     (n1->ndx + 2 == n2->ndx || n1->ndx - 1 == n2->ndx)) ||
    (n1->strand == -1 && n2->strand == -1 &&
     (n1->ndx + 2 == n2->ndx || n1->ndx - 1 == n2->ndx))) {
    if(n1->strand == 1 && n2->rscore < 0) rval -= n2->rscore;
    if(n1->strand == -1 && n1->rscore < 0) rval -= n1->rscore;
    if(n1->strand == 1 && n2->uscore < 0) rval -= n2->uscore;
    if(n1->strand == -1 && n1->uscore < 0) rval -= n1->uscore;
  }
  dist = abs(n1->ndx-n2->ndx);
  if(n1->strand == 1 && n2->strand == 1 && n1->ndx+2 >= n2->ndx) ovlp = 1;
  else if(n1->strand == -1 && n2->strand == -1 && n1->ndx >= n2->ndx+2)
    ovlp = 1;
  if(dist > 3*OPER_DIST || n1->strand != n2->strand) {
    rval -= 0.15 * start_weight;
  }
  else if((dist <= OPER_DIST && ovlp == 0) || dist < 0.25*OPER_DIST) {
    rval += (2.0 - (double)(dist)/OPER_DIST) * 0.15 * start_weight;
  }
  return rval;
}

/* Return the minimum of two numbers */

double dmin(double x, double y) {
  if(x < y) return x;
  return y;
}

/* Return the maximum of two numbers */

double dmax(double x, double y) {
  if(x > y) return x;
  return y;
}

/* Sorting routine for nodes */

int compare_nodes(const void *v1, const void *v2) {
  struct _node *n1, *n2;
  n1 = (struct _node *)v1;
  n2 = (struct _node *)v2;
  if(n1->ndx < n2->ndx) return -1;
  if(n1->ndx > n2->ndx) return 1;
  if(n1->strand > n2->strand) return -1;
  if(n1->strand < n2->strand) return 1;
  return 0;
}

/* Sorts all nodes by common stop */

int stopcmp_nodes(const void *v1, const void *v2) {
  struct _node *n1, *n2;
  n1 = (struct _node *)v1;
  n2 = (struct _node *)v2;
  if(n1->stop_val < n2->stop_val) return -1;
  if(n1->stop_val > n2->stop_val) return 1;
  if(n1->strand > n2->strand) return -1;
  if(n1->strand < n2->strand) return 1;
  if(n1->ndx < n2->ndx) return -1;
  if(n1->ndx > n2->ndx) return 1;
  return 0;
}
