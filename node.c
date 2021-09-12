/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

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
  run off the edge, in which case they only have to be 50bp.
*******************************************************************************/

int add_nodes(unsigned char *seq, unsigned char *rseq, int slen, struct _node
              *nodes, int closed, mask *mlist, int nm, struct _training
              *tinf) {
  int i, nn = 0, last[3], saw_start[3], min_dist[3];
  int slmod = 0;

  /* Forward strand nodes */
  slmod = slen%3;
  for(i = 0; i < 3; i++) {
    last[(i+slmod)%3] = slen+i; 
    saw_start[i%3] = 0;
    min_dist[i%3] = MIN_EDGE_GENE;
    if(closed == 0) while(last[(i+slmod)%3]+2 > slen-1) last[(i+slmod)%3]-=3;
  }
  for(i = slen-3; i >= 0; i--) {
    if(is_stop(seq, i, tinf)==1) {
      if(saw_start[i%3] == 1) {
        if(is_stop(seq, last[i%3], tinf) == 0) nodes[nn].edge = 1;
        nodes[nn].ndx = last[i%3]; 
        nodes[nn].type = STOP;
        nodes[nn].strand = 1; 
        nodes[nn++].stop_val = i;
      }
      min_dist[i%3] = MIN_GENE;
      last[i%3]=i; 
      saw_start[i%3] = 0;
      continue;
    }
    if(last[i%3] >= slen) continue;
     
    if(is_start(seq, i, tinf) == 1 && is_atg(seq, i)==1 && ((last[i%3]-i+3)
            >= min_dist[i%3]) && cross_mask(i, last[i%3], mlist, nm) == 0) {
      nodes[nn].ndx = i; 
      nodes[nn].type = ATG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = last[i%3]; 
      nodes[nn++].strand = 1;
    }
    else if(is_start(seq, i, tinf) == 1 && is_gtg(seq, i)==1 && ((last[i%3]-i+3)
            >= min_dist[i%3]) && cross_mask(i, last[i%3], mlist, nm) == 0) {
      nodes[nn].ndx = i; 
      nodes[nn].type = GTG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = last[i%3]; 
      nodes[nn++].strand = 1;
    }
    else if(is_start(seq, i, tinf) == 1 && is_ttg(seq, i)==1 && ((last[i%3]-i+3)
            >= min_dist[i%3]) && cross_mask(i, last[i%3], mlist, nm) == 0) {
      nodes[nn].ndx = i; 
      nodes[nn].type = TTG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = last[i%3]; 
      nodes[nn++].strand = 1;
    }
    else if(i <= 2 && closed == 0 && ((last[i%3]-i) > MIN_EDGE_GENE) &&
            cross_mask(i, last[i%3], mlist, nm) == 0) {
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
      if(is_stop(seq, last[i%3], tinf) == 0) nodes[nn].edge = 1;
      nodes[nn].ndx = last[i%3]; 
      nodes[nn].type = STOP;
      nodes[nn].strand = 1; 
      nodes[nn++].stop_val = i-6;
    }
  }

  /* Reverse strand nodes */
  for(i = 0; i < 3; i++) {
    last[(i+slmod)%3] = slen+i; 
    saw_start[i%3] = 0;
    min_dist[i%3] = MIN_EDGE_GENE;
    if(closed == 0) while(last[(i+slmod)%3]+2 > slen-1) last[(i+slmod)%3]-=3;
  }
  for(i = slen-3; i >= 0; i--) {
    if(is_stop(rseq, i, tinf)==1) {
      if(saw_start[i%3] == 1) {
        if(is_stop(rseq, last[i%3], tinf) == 0) nodes[nn].edge = 1;
        nodes[nn].ndx = slen-last[i%3]-1; 
        nodes[nn].type = STOP;
        nodes[nn].strand = -1; 
        nodes[nn++].stop_val = slen-i-1;
      }
      min_dist[i%3] = MIN_GENE;
      last[i%3]=i; 
      saw_start[i%3] = 0;
      continue;
    }
    if(last[i%3] >= slen) continue;

    if(is_start(rseq, i, tinf) == 1 && is_atg(rseq, i)==1 && ((last[i%3]-i+3)
       >= min_dist[i%3]) && cross_mask(slen-last[i%3]-1, slen-i-1, mlist, nm) ==
       0) {
      nodes[nn].ndx = slen - i - 1; 
      nodes[nn].type = ATG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = slen-last[i%3]-1; 
      nodes[nn++].strand = -1;
    }
    else if(is_start(rseq, i, tinf) == 1 && is_gtg(rseq, i)==1 && 
            ((last[i%3]-i+3) >= min_dist[i%3]) && cross_mask(slen-last[i%3]-1,
            slen-i-1, mlist, nm) == 0) {
      nodes[nn].ndx = slen - i - 1; 
      nodes[nn].type = GTG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = slen-last[i%3]-1; 
      nodes[nn++].strand = -1;
    }
    else if(is_start(rseq, i, tinf) == 1 && is_ttg(rseq, i)==1 && 
            ((last[i%3]-i+3) >= min_dist[i%3]) && cross_mask(slen-last[i%3]-1,
            slen-i-1, mlist, nm) == 0) {
      nodes[nn].ndx = slen - i - 1; 
      nodes[nn].type = TTG; 
      saw_start[i%3] = 1;
      nodes[nn].stop_val = slen-last[i%3]-1; 
      nodes[nn++].strand = -1;
    }
    else if(i <= 2 && closed == 0 && ((last[i%3]-i) > MIN_EDGE_GENE) &&
            cross_mask(slen-last[i%3]-1, slen-i-1, mlist, nm) == 0) {
      nodes[nn].ndx = slen - i - 1; 
      nodes[nn].type = ATG; 
      saw_start[i%3] = 1;
      nodes[nn].edge = 1; 
      nodes[nn].stop_val = slen-last[i%3]-1;
      nodes[nn++].strand = -1;
    }
  }
  for(i = 0; i < 3; i++) {
    if(saw_start[i%3] == 1) {
      if(is_stop(rseq, last[i%3], tinf) == 0) nodes[nn].edge = 1;
      nodes[nn].ndx = slen - last[i%3] - 1; 
      nodes[nn].type = STOP;
      nodes[nn].strand = -1; 
      nodes[nn++].stop_val = slen-i+5;
    }
  }
  return nn;
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

void record_overlapping_starts(struct _node *nod, int nn, struct _training
                               *tinf, int flag) {
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
                  intergenic_mod(&nod[i], &nod[j], tinf) > max_sc)) {
            nod[i].star_ptr[(nod[j].ndx)%3] = j;
            max_sc = nod[j].cscore + nod[j].sscore +
                     intergenic_mod(&nod[i], &nod[j], tinf);
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
                  intergenic_mod(&nod[j], &nod[i], tinf) > max_sc)) {
            nod[i].star_ptr[(nod[j].ndx)%3] = j;
            max_sc = nod[j].cscore + nod[j].sscore +
                     intergenic_mod(&nod[j], &nod[i], tinf);
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

void record_gc_bias(int *gc, struct _node *nod, int nn, struct _training
                    *tinf) {
  int i, j, ctr[3][3], last[3], frmod, fr, mfr, len;
  double tot = 0.0;

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

  for(i = 0; i < 3; i++) tinf->bias[i] = 0.0;
  for(i = 0; i < nn; i++) {
    if(nod[i].type != STOP) {
      len = abs(nod[i].stop_val-nod[i].ndx)+1;
      tinf->bias[nod[i].gc_bias]+= (nod[i].gc_score[nod[i].gc_bias]*len)/1000.0;
    }
  }
  tot = tinf->bias[0] + tinf->bias[1] + tinf->bias[2];
  for(i = 0; i < 3; i++) tinf->bias[i] *= (3.0/tot);
}

/*******************************************************************************
  Simple routine that calculates the dicodon frequency in genes and in the
  background, and then stores the log likelihood of each 6-mer relative to the
  background.
*******************************************************************************/

void calc_dicodon_gene(struct _training *tinf, unsigned char *seq, unsigned
                       char *rseq, int slen, struct _node *nod, int dbeg) {
  int i, path, counts[4096], glob = 0;
  int left, right, in_gene;
  double prob[4096], bg[4096];

  for(i = 0; i < 4096; i++) { counts[i] = 0; prob[i] = 0.0; bg[i] = 0.0; }
  left = -1; right = -1;
  calc_mer_bg(6, seq, rseq, slen, bg);
  path = dbeg; in_gene = 0;
  while(path != -1) {
    if(nod[path].strand == -1 && nod[path].type != STOP) {
      in_gene = -1;
      left = slen-nod[path].ndx-1;
    }
    if(nod[path].strand == 1 && nod[path].type == STOP) {
      in_gene = 1;
      right = nod[path].ndx+2;
    }
    if(in_gene == -1 && nod[path].strand == -1 && nod[path].type == STOP) {
      right = slen-nod[path].ndx+1;
      for(i = left; i < right-5; i+=3) {
        counts[mer_ndx(6, rseq, i)]++;
        glob++;
      }
      in_gene = 0;
    }
    if(in_gene == 1 && nod[path].strand == 1 && nod[path].type != STOP) {
      left = nod[path].ndx;
      for(i = left; i < right-5; i+=3) { counts[mer_ndx(6, seq, i)]++; glob++; }
      in_gene = 0;
    }
    path = nod[path].traceb;
  }
  for(i = 0; i < 4096; i++) {
    prob[i] = (counts[i]*1.0)/(glob*1.0);
    if(prob[i] == 0 && bg[i] != 0) tinf->gene_dc[i] = -5.0;
    else if(bg[i] == 0) tinf->gene_dc[i] = 0.0;
    else tinf->gene_dc[i] = log(prob[i]/bg[i]);
    if(tinf->gene_dc[i] > 5.0) tinf->gene_dc[i] = 5.0;
    if(tinf->gene_dc[i] < -5.0) tinf->gene_dc[i] = -5.0;
  }
}

/*******************************************************************************
  Scoring function for all the start nodes.  This score has two factors:  (1)
  Coding, which is a composite of coding score and length, and (2) Start
  score, which is a composite of RBS score and ATG/TTG/GTG.
*******************************************************************************/

void score_nodes(unsigned char *seq, unsigned char *rseq, int slen,
                 struct _node *nod, int nn, struct _training *tinf,
                 int closed, int is_meta) {
  int i, j;
  double negf, posf, rbs1, rbs2, sd_score, edge_gene, min_meta_len;

  /* Step 1: Calculate raw coding potential for every start-stop pair. */
  calc_orf_gc(seq, rseq, slen, nod, nn, tinf);
  raw_coding_score(seq, rseq, slen, nod, nn, tinf);

  /* Step 2: Calculate raw RBS Scores for every start node. */
  if(tinf->uses_sd == 1) rbs_score(seq, rseq, slen, nod, nn, tinf);
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
       tinf) == 0) || (nod[i].strand == -1 && is_stop(rseq, slen-1-
       nod[i].stop_val, tinf) == 0)) edge_gene++;

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
    /* Coding Penalization in Metagenomic Fragments:  Internal    */
    /* genes must have a score of 5.0 and be >= 120bp.  High GC   */
    /* genes are also penalized.                                  */
    /**************************************************************/
    if(is_meta == 1 && slen < 3000 && edge_gene == 0 && 
       (nod[i].cscore < 5.0 || abs(nod[i].ndx-nod[i].stop_val) < 120)) {
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
        if(is_meta == 0 || slen > 1500) nod[i].sscore -= tinf->st_wt;
        else nod[i].sscore -= (10.31 - 0.004*slen);
      }
      else if(is_meta == 1 && slen < 3000 && nod[i].edge == 1) {
        min_meta_len = sqrt(slen)*5.0;
        if(abs(nod[i].ndx-nod[i].stop_val) >= min_meta_len) {
          if(nod[i].cscore >= 0) nod[i].cscore = -1.0;
          nod[i].sscore = 0.0; 
          nod[i].uscore = 0.0; 
        }
      }
      else nod[i].sscore -= 0.5;
    }
    else if(nod[i].cscore < 5.0 && is_meta == 1 && abs(nod[i].ndx-
            nod[i].stop_val) < 120 && nod[i].sscore < 0.0)
      nod[i].sscore -= tinf->st_wt; 
  }
}

/* Calculate the GC Content for each start-stop pair */
void calc_orf_gc(unsigned char *seq, unsigned char *rseq, int slen, struct
                 _node *nod, int nn, struct _training *tinf) {
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
                      _node *nod, int nn, struct _training *tinf) {
  int i, j, last[3], fr;
  double score[3], lfac, no_stop, gsize = 0.0;

  if(tinf->trans_table != 11) { /* TGA or TAG is not a stop */
    no_stop = ((1-tinf->gc)*(1-tinf->gc)*tinf->gc)/8.0;
    no_stop += ((1-tinf->gc)*(1-tinf->gc)*(1-tinf->gc))/8.0;
    no_stop = (1 - no_stop);
  }
  else {
    no_stop = ((1-tinf->gc)*(1-tinf->gc)*tinf->gc)/4.0;
    no_stop += ((1-tinf->gc)*(1-tinf->gc)*(1-tinf->gc))/8.0;
    no_stop = (1 - no_stop);
  }

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
        score[fr] += tinf->gene_dc[mer_ndx(6, seq, j)];
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
        score[fr] += tinf->gene_dc[mer_ndx(6, rseq, slen-j-1)];
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
  Examines the results of the SD motif search to determine if this organism
  uses an SD motif or not.  Some motif of 3-6bp has to be good or we set
  uses_sd to 0, which will cause Prodigal to run the non-SD motif finder for
  starts.
*******************************************************************************/
void determine_sd_usage(struct _training *tinf) {
  tinf->uses_sd = 1;
  if(tinf->rbs_wt[0] >= 0.0) tinf->uses_sd = 0;
  if(tinf->rbs_wt[16] < 1.0 && tinf->rbs_wt[13] < 1.0 && tinf->rbs_wt[15] < 1.0
     && (tinf->rbs_wt[0] >= -0.5 || (tinf->rbs_wt[22] < 2.0 && tinf->rbs_wt[24]
     < 2.0 && tinf->rbs_wt[27] < 2.0)))
    tinf->uses_sd = 0;
}

/*******************************************************************************
  RBS Scoring Function: Calculate the RBS motif and then multiply it by the
  appropriate weight for that motif (determined in the start training
  function).
*******************************************************************************/
void rbs_score(unsigned char *seq, unsigned char *rseq, int slen, struct _node
               *nod, int nn, struct _training *tinf) {
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
        cur_sc[0] = shine_dalgarno_exact(seq, j, nod[i].ndx, tinf->rbs_wt);
        cur_sc[1] = shine_dalgarno_mm(seq, j, nod[i].ndx, tinf->rbs_wt);
        if(cur_sc[0] > nod[i].rbs[0]) nod[i].rbs[0] = cur_sc[0];
        if(cur_sc[1] > nod[i].rbs[1]) nod[i].rbs[1] = cur_sc[1];
      }
    }
    else if(nod[i].strand == -1) {
      for(j = slen - nod[i].ndx - 21; j <= slen - nod[i].ndx - 7; j++) {
        if(j > slen-1) continue;
        cur_sc[0] = shine_dalgarno_exact(rseq, j, slen-1-nod[i].ndx,
                                         tinf->rbs_wt);
        cur_sc[1] = shine_dalgarno_mm(rseq, j, slen-1-nod[i].ndx,
                                      tinf->rbs_wt);
        if(cur_sc[0] > nod[i].rbs[0]) nod[i].rbs[0] = cur_sc[0];
        if(cur_sc[1] > nod[i].rbs[1]) nod[i].rbs[1] = cur_sc[1];
      }
    }
  }
}

/*******************************************************************************
  Iterative Algorithm to train starts.  It begins with all the highest coding
  starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts moving
  starts around attempting to match these discoveries.  This start trainer is
  for Shine-Dalgarno motifs only.
*******************************************************************************/
void train_starts_sd(unsigned char *seq, unsigned char *rseq, int slen,
                  struct _node *nod, int nn, struct _training *tinf) {
  int i, j, fr, rbs[3], type[3], bndx[3], max_rb;
  double sum, wt, rbg[28], rreal[28], best[3], sthresh = 35.0;
  double tbg[3], treal[3];

  wt = tinf->st_wt;
  for(j = 0; j < 3; j++) tinf->type_wt[j] = 0.0;
  for(j = 0; j < 28; j++) tinf->rbs_wt[j] = 0.0;
  for(i = 0; i < 32; i++) for(j = 0; j < 4; j++) tinf->ups_comp[i][j] = 0.0;

  /* Build the background of random types */
  for(i = 0; i < 3; i++) tbg[i] = 0.0;
  for(i = 0; i < nn; i++) {
    if(nod[i].type == STOP) continue;
    tbg[nod[i].type] += 1.0;
  }
  sum = 0.0;
  for(i = 0; i < 3; i++) sum += tbg[i];
  for(i = 0; i < 3; i++) tbg[i] /= sum;

  /* Iterate 10 times through the list of nodes                         */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few     */
  /* extra to be safe)                                                  */
  for(i = 0; i < 10; i++) {

    /* Recalculate the RBS motif background */
    for(j = 0; j < 28; j++) rbg[j] = 0.0;
    for(j = 0; j < nn; j++) {
      if(nod[j].type == STOP || nod[j].edge == 1) continue;
      if(tinf->rbs_wt[nod[j].rbs[0]] > tinf->rbs_wt[nod[j].rbs[1]]+1.0 ||
         nod[j].rbs[1] == 0)
        max_rb = nod[j].rbs[0];
      else if(tinf->rbs_wt[nod[j].rbs[0]] < tinf->rbs_wt[nod[j].rbs[1]]-1.0 ||
        nod[j].rbs[0] == 0)
        max_rb = nod[j].rbs[1];
      else max_rb = (int)dmax(nod[j].rbs[0], nod[j].rbs[1]);
      rbg[max_rb] += 1.0;
    }
    sum = 0.0;
    for(j = 0; j < 28; j++) sum += rbg[j];
    for(j = 0; j < 28; j++) rbg[j] /= sum;

    for(j = 0; j < 28; j++) rreal[j] = 0.0;
    for(j = 0; j < 3; j++) treal[j] = 0.0;

    /* Forward strand pass */
    for(j = 0; j < 3; j++) {
      best[j] = 0.0; bndx[j] = -1; rbs[j] = 0; type[j] = 0;
    }
    for(j = 0; j < nn; j++) {
      if(nod[j].type != STOP && nod[j].edge == 1) continue;
      fr = (nod[j].ndx)%3;
      if(nod[j].type == STOP && nod[j].strand == 1) {
        if(best[fr] >= sthresh && nod[bndx[fr]].ndx%3 == fr) {
          rreal[rbs[fr]] += 1.0;
          treal[type[fr]] += 1.0;
          if(i == 9) count_upstream_composition(seq, slen, 1,
                      nod[bndx[fr]].ndx, tinf);          
        }
        best[fr] = 0.0; bndx[fr] = -1; rbs[fr] = 0; type[fr] = 0;
      }
      else if(nod[j].strand == 1) {
        if(tinf->rbs_wt[nod[j].rbs[0]] > tinf->rbs_wt[nod[j].rbs[1]]+1.0 ||
           nod[j].rbs[1] == 0)
          max_rb = nod[j].rbs[0];
        else if(tinf->rbs_wt[nod[j].rbs[0]] < tinf->rbs_wt[nod[j].rbs[1]]-1.0 ||
                nod[j].rbs[0] == 0)
          max_rb = nod[j].rbs[1];
        else max_rb = (int)dmax(nod[j].rbs[0], nod[j].rbs[1]);
        if(nod[j].cscore + wt*tinf->rbs_wt[max_rb] +
           wt*tinf->type_wt[nod[j].type] >= best[fr]) {
          best[fr] = nod[j].cscore + wt*tinf->rbs_wt[max_rb];
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bndx[fr] = j;
          type[fr] = nod[j].type;
          rbs[fr] = max_rb;
        }
      }
    }

    /* Reverse strand pass */
    for(j = 0; j < 3; j++) {
      best[j] = 0.0; bndx[j] = -1; rbs[j] = 0; type[j] = 0;
    }
    for(j = nn-1; j >= 0; j--) {
      if(nod[j].type != STOP && nod[j].edge == 1) continue;
      fr = (nod[j].ndx)%3;
      if(nod[j].type == STOP && nod[j].strand == -1) {
        if(best[fr] >= sthresh && nod[bndx[fr]].ndx%3 == fr) {
          rreal[rbs[fr]] += 1.0;
          treal[type[fr]] += 1.0;
          if(i == 9) count_upstream_composition(rseq, slen, -1,
                      nod[bndx[fr]].ndx, tinf);          
        }
        best[fr] = 0.0; bndx[fr] = -1; rbs[fr] = 0; type[fr] = 0;
      }
      else if(nod[j].strand == -1) {
        if(tinf->rbs_wt[nod[j].rbs[0]] > tinf->rbs_wt[nod[j].rbs[1]]+1.0 ||
           nod[j].rbs[1] == 0)
          max_rb = nod[j].rbs[0];
        else if(tinf->rbs_wt[nod[j].rbs[0]] < tinf->rbs_wt[nod[j].rbs[1]]-1.0 ||
           nod[j].rbs[0] == 0)
          max_rb = nod[j].rbs[1];
        else max_rb = (int)dmax(nod[j].rbs[0], nod[j].rbs[1]);
        if(nod[j].cscore + wt*tinf->rbs_wt[max_rb] +
           wt*tinf->type_wt[nod[j].type] >= best[fr]) {
          best[fr] = nod[j].cscore + wt*tinf->rbs_wt[max_rb];
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bndx[fr] = j;
          type[fr] = nod[j].type;
          rbs[fr] = max_rb;
        }
      }
    }

    sum = 0.0;
    for(j = 0; j < 28; j++) sum += rreal[j];
    if(sum == 0.0) for(j = 0; j < 28; j++) tinf->rbs_wt[j] = 0.0;
    else {
      for(j = 0; j < 28; j++) {
        rreal[j] /= sum;
        if(rbg[j] != 0) tinf->rbs_wt[j] = log(rreal[j]/rbg[j]);
        else tinf->rbs_wt[j] = -4.0;
        if(tinf->rbs_wt[j] > 4.0) tinf->rbs_wt[j] = 4.0;
        if(tinf->rbs_wt[j] < -4.0) tinf->rbs_wt[j] = -4.0;
      }
    }
    sum = 0.0;
    for(j = 0; j < 3; j++) sum += treal[j];
    if(sum == 0.0) for(j = 0; j < 3; j++) tinf->type_wt[j] = 0.0;
    else {
      for(j = 0; j < 3; j++) {
        treal[j] /= sum;
        if(tbg[j] != 0) tinf->type_wt[j] = log(treal[j]/tbg[j]);
        else tinf->type_wt[j] = -4.0;
        if(tinf->type_wt[j] > 4.0) tinf->type_wt[j] = 4.0;
        if(tinf->type_wt[j] < -4.0) tinf->type_wt[j] = -4.0;
      }
    }
    if(sum <= (double)nn/2000.0) sthresh /= 2.0;
  }

/* Convert upstream base composition to a log score */
for(i = 0; i < 32; i++) {
  sum = 0.0;
  for(j = 0; j < 4; j++) sum += tinf->ups_comp[i][j];
  if(sum == 0.0) for(j = 0; j < 4; j++) tinf->ups_comp[i][j] = 0.0;
  else {
    for(j = 0; j < 4; j++) {
      tinf->ups_comp[i][j] /= sum;
      if(tinf->gc > 0.1 && tinf->gc < 0.9) {
        if(j == 0 || j == 3)
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/(1.0-tinf->gc));
        else
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/tinf->gc);
      }
      else if(tinf->gc <= 0.1) {
        if(j == 0 || j == 3) 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
        else 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
      }
      else {
        if(j == 0 || j == 3) 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
        else 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
      }
      if(tinf->ups_comp[i][j] > 4.0) tinf->ups_comp[i][j] = 4.0;
      if(tinf->ups_comp[i][j] < -4.0) tinf->ups_comp[i][j] = -4.0;
    }     
  }
}

/* Start training info: kept this in since it's useful information */
/* to print sometimes. */
/* fprintf(stderr, "\nLOG WTS\n");
for(i = 0; i < 3; i++) fprintf(stderr, "%f ", tinf->type_wt[i]);
fprintf(stderr, "\n");
for(i = 0; i < 28; i++) fprintf(stderr, "%f ", tinf->rbs_wt[i]);
fprintf(stderr, "\n\nSTART DIST: ");
for(i = 0; i < 3; i++) fprintf(stderr, "%f ", treal[i]);
fprintf(stderr, "\n");
sum = 0.0;
for(i = 0; i < 28; i++) { fprintf(stderr, "%f ", rreal[i]); sum+= rreal[i]; }
fprintf(stderr, "sum is %f\n", sum);
fprintf(stderr, "\n\nUPS COMP: ");
for(i = 0; i < 32; i++) { fprintf(stderr, "%d", i); for(j = 0; j < 4; j++) { fprintf(stderr, "\t%.2f", tinf->ups_comp[i][j]); } fprintf(stderr, "\n"); } 
exit(0); */
}

/*******************************************************************************
  Iterative Algorithm to train starts.  It begins with all the highest coding
  starts in the model, scans for RBS/ATG-GTG-TTG usage, then starts moving
  starts around attempting to match these discoveries.  Unlike the SD
  algorithm, it allows for any popular motif to be discovered.
*******************************************************************************/
void train_starts_nonsd(unsigned char *seq, unsigned char *rseq, int slen,
                  struct _node *nod, int nn, struct _training *tinf) {
  int i, j, k, l, fr, bndx[3], mgood[4][4][4096], stage;
  double sum, ngenes, wt = tinf->st_wt, best[3], sthresh = 35.0;
  double tbg[3], treal[3];
  double mbg[4][4][4096], mreal[4][4][4096], zbg, zreal;

  for(i = 0; i < 32; i++) for(j = 0; j < 4; j++) tinf->ups_comp[i][j] = 0.0;

  /* Build the background of random types */
  for(i = 0; i < 3; i++) tinf->type_wt[i] = 0.0;
  for(i = 0; i < 3; i++) tbg[i] = 0.0;
  for(i = 0; i < nn; i++) {
    if(nod[i].type == STOP) continue;
    tbg[nod[i].type] += 1.0;
  }
  sum = 0.0;
  for(i = 0; i < 3; i++) sum += tbg[i];
  for(i = 0; i < 3; i++) tbg[i] /= sum;

  /* Iterate 20 times through the list of nodes                         */
  /* Converge upon optimal weights for ATG vs GTG vs TTG and RBS motifs */
  /* (convergence typically takes 4-5 iterations, but we run a few      */
  /* extra to be safe)                                                  */
  for(i = 0; i < 20; i++) {

    /* Determine which stage of motif finding we're in */
    if(i < 4) stage = 0;
    else if(i < 12) stage = 1;
    else stage = 2;

    /* Recalculate the upstream motif background and set 'real' counts to 0 */
    for(j = 0; j < 4; j++) for(k = 0; k < 4; k++) for(l = 0; l < 4096; l++)
      mbg[j][k][l] = 0.0;
    zbg = 0.0;
    for(j = 0; j < nn; j++) {
      if(nod[j].type == STOP || nod[j].edge == 1) continue;
      find_best_upstream_motif(tinf, seq, rseq, slen, &nod[j], stage);
      update_motif_counts(mbg, &zbg, seq, rseq, slen, &(nod[j]), stage);
    }
    sum = 0.0;
    for(j = 0; j < 4; j++) for(k = 0; k < 4; k++) for(l = 0; l < 4096; l++)
      sum += mbg[j][k][l];
    sum += zbg;
    for(j = 0; j < 4; j++) for(k = 0; k < 4; k++) for(l = 0; l < 4096; l++)
      mbg[j][k][l] /= sum;
    zbg /= sum;

    /* Reset counts of 'real' motifs/types to 0 */
    for(j = 0; j < 4; j++) for(k = 0; k < 4; k++) for(l = 0; l < 4096; l++)
      mreal[j][k][l] = 0.0;
    zreal = 0.0;
    for(j = 0; j < 3; j++) treal[j] = 0.0;
    ngenes = 0.0;

    /* Forward strand pass */
    for(j = 0; j < 3; j++) { best[j] = 0.0; bndx[j] = -1; }
    for(j = 0; j < nn; j++) {
      if(nod[j].type != STOP && nod[j].edge == 1) continue;
      fr = (nod[j].ndx)%3;
      if(nod[j].type == STOP && nod[j].strand == 1) {
        if(best[fr] >= sthresh) {
          ngenes += 1.0;
          treal[nod[bndx[fr]].type] += 1.0;
          update_motif_counts(mreal, &zreal, seq, rseq, slen, &(nod[bndx[fr]]),
                              stage);
          if(i == 19) count_upstream_composition(seq, slen, 1,
                      nod[bndx[fr]].ndx, tinf);          
        }
        best[fr] = 0.0; bndx[fr] = -1;
      }
      else if(nod[j].strand == 1) {
        if(nod[j].cscore + wt*nod[j].mot.score + wt*tinf->type_wt[nod[j].type]
           >= best[fr]) {
          best[fr] = nod[j].cscore + wt*nod[j].mot.score;
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bndx[fr] = j;
        }
      }
    }

    /* Reverse strand pass */
    for(j = 0; j < 3; j++) { best[j] = 0.0; bndx[j] = -1; }
    for(j = nn-1; j >= 0; j--) {
      if(nod[j].type != STOP && nod[j].edge == 1) continue;
      fr = (nod[j].ndx)%3;
      if(nod[j].type == STOP && nod[j].strand == -1) {
        if(best[fr] >= sthresh) {
          ngenes += 1.0;
          treal[nod[bndx[fr]].type] += 1.0;
          update_motif_counts(mreal, &zreal, seq, rseq, slen, &(nod[bndx[fr]]),
                              stage);
          if(i == 19) count_upstream_composition(rseq, slen, -1,
                      nod[bndx[fr]].ndx, tinf);          
        }
        best[fr] = 0.0; bndx[fr] = -1;
      }
      else if(nod[j].strand == -1) {
        if(nod[j].cscore + wt*nod[j].mot.score + wt*tinf->type_wt[nod[j].type]
           >= best[fr]) {
          best[fr] = nod[j].cscore + wt*nod[j].mot.score;
          best[fr] += wt*tinf->type_wt[nod[j].type];
          bndx[fr] = j;
        }
      }
    }

    /* Update the log likelihood weights for type and RBS motifs */
    if(stage < 2) build_coverage_map(mreal, mgood, ngenes, stage);
    sum = 0.0;
    for(j = 0; j < 4; j++) for(k = 0; k < 4; k++) for(l = 0; l < 4096; l++)
      sum += mreal[j][k][l];
    sum += zreal;
    if(sum == 0.0) {
      for(j = 0; j < 4; j++) for(k = 0; k < 4; k++) for(l = 0; l < 4096; l++)
        tinf->mot_wt[j][k][l] = 0.0;
      tinf->no_mot = 0.0;
    }
    else {
      for(j = 0; j < 4; j++) for(k = 0; k < 4; k++)
      for(l = 0; l < 4096; l++) {{{
        if(mgood[j][k][l] == 0) {
          zreal += mreal[j][k][l];
          zbg += mreal[j][k][l];
          mreal[j][k][l] = 0.0;
          mbg[j][k][l] = 0.0;
        }
        mreal[j][k][l] /= sum;
        if(mbg[j][k][l] != 0)
          tinf->mot_wt[j][k][l] = log(mreal[j][k][l]/mbg[j][k][l]);
        else tinf->mot_wt[j][k][l] = -4.0;
        if(tinf->mot_wt[j][k][l] > 4.0) tinf->mot_wt[j][k][l] = 4.0;
        if(tinf->mot_wt[j][k][l] < -4.0) tinf->mot_wt[j][k][l] = -4.0;
      }}}
    }
    zreal /= sum;
    if(zbg != 0) tinf->no_mot = log(zreal/zbg);
    else tinf->no_mot = -4.0;
    if(tinf->no_mot > 4.0) tinf->no_mot = 4.0;
    if(tinf->no_mot < -4.0) tinf->no_mot = -4.0;
    sum = 0.0;
    for(j = 0; j < 3; j++) sum += treal[j];
    if(sum == 0.0) for(j = 0; j < 3; j++) tinf->type_wt[j] = 0.0;
    else {
      for(j = 0; j < 3; j++) {
        treal[j] /= sum;
        if(tbg[j] != 0) tinf->type_wt[j] = log(treal[j]/tbg[j]);
        else tinf->type_wt[j] = -4.0;
        if(tinf->type_wt[j] > 4.0) tinf->type_wt[j] = 4.0;
        if(tinf->type_wt[j] < -4.0) tinf->type_wt[j] = -4.0;
      }
    }
    if(sum <= (double)nn/2000.0) sthresh /= 2.0;
  }

/* Convert upstream base composition to a log score */
for(i = 0; i < 32; i++) {
  sum = 0.0;
  for(j = 0; j < 4; j++) sum += tinf->ups_comp[i][j];
  if(sum == 0.0) for(j = 0; j < 4; j++) tinf->ups_comp[i][j] = 0.0;
  else {
    for(j = 0; j < 4; j++) {
      tinf->ups_comp[i][j] /= sum;
      if(tinf->gc > 0.1 && tinf->gc < 0.9) {
        if(j == 0 || j == 3)
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/(1.0-tinf->gc));
        else
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/tinf->gc);
      }
      else if(tinf->gc <= 0.1) {
        if(j == 0 || j == 3) 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
        else 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
      }
      else {
        if(j == 0 || j == 3) 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.10);
        else 
          tinf->ups_comp[i][j] = log(tinf->ups_comp[i][j]*2.0/0.90);
      }
      if(tinf->ups_comp[i][j] > 4.0) tinf->ups_comp[i][j] = 4.0;
      if(tinf->ups_comp[i][j] < -4.0) tinf->ups_comp[i][j] = -4.0;
    }     
  }
}

/* fprintf(stderr, "\n\nUPS COMP: ");
for(i = 0; i < 32; i++) { fprintf(stderr, "%d", i); for(j = 0; j < 4; j++) { fprintf(stderr, "\t%.2f", tinf->ups_comp[i][j]); } fprintf(stderr, "\n"); }  exit(0);
*/
/* motif dump, keeping it in since it can be useful
for(j = 0; j < 4096; j++) {
  for(k = 0; k < 4; k++) {
    if(k == 0 && j >= 64) continue;
    if(k == 1 && j >= 256) continue;
    if(k == 2 && j >= 1025) continue;
    mer_text(qt, k+3, j);
    printf("%s\t", qt);
    for(l = 0; l < 4; l++) {
      printf("Dist %d Counts %.2f BG %.2f Val %.2f\t", l, mreal[k][l][j],
             mbg[k][l][j], tinf->mot_wt[k][l][j]);
    }
    printf("\n");
  }
}
exit(0);
*/

}

/*******************************************************************************
  For a given start, record the base composition of the upstream region at
  positions -1 and -2 and -15 to -44.  This will be used to supplement the
  SD (or other) motif finder with additional information.
*******************************************************************************/
void count_upstream_composition(unsigned char *seq, int slen, int strand, 
                                int pos, struct _training *tinf) {
  int i, start, count = 0;
  if(strand == 1) start = pos; 
  else start = slen-1-pos; 

  for(i = 1; i < 45; i++) {
    if(i > 2 && i < 15) continue;
    if(start-i >= 0) tinf->ups_comp[count][mer_ndx(1, seq, start-i)]++;
    count++;
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
  Update the motif counts from a putative "real" start.  This is done in three
  stages.  In stage 0, all motifs sizes 3-6bp in the region with spacer 3-15bp
  are counted.  In stage 1, only the best motif and all its subsets are
  counted (e.g. for AGGAG, we would count AGGAG, AGGA, GGAG, AGG, GGA, and
  GAG).  In stage 2, only the best single motif is counted.
*******************************************************************************/
void update_motif_counts(double mcnt[4][4][4096], double *zero, unsigned char
                         *seq, unsigned char *rseq, int slen, struct _node *nod,
                         int stage) {
  int i, j, k, start, spacendx;
  unsigned char *wseq;
  struct _motif *mot = &(nod->mot);

  if(nod->type == STOP || nod->edge == 1) return;
  if(mot->len == 0) { *zero += 1.0; return; }

  if(nod->strand == 1) { wseq = seq; start = nod->ndx; }
  else { wseq = rseq; start = slen-1-nod->ndx; }

  /* Stage 0:  Count all motifs.  If a motif is detected, */
  /* it is counted for every distance in stage 0.  This   */
  /* is done to make sure off-distance good motifs are    */
  /* recognized.                                          */
  if(stage == 0) {
    for(i = 3; i >= 0; i--) {
      for(j = start-18-i; j <= start-6-i; j++) {
        if(j < 0) continue;
        if(j <= start-16-i) spacendx = 3;
        else if(j <= start-14-i) spacendx = 2;
        else if(j >= start-7-i) spacendx = 1;
        else spacendx = 0;
        for(k = 0; k < 4; k++) mcnt[i][k][mer_ndx(i+3, wseq, j)] += 1.0;
      }
    }
  }
  /* Stage 1:  Count only the best motif, but also count  */
  /* all its sub-motifs.                                  */
  else if(stage == 1) {
    mcnt[mot->len-3][mot->spacendx][mot->ndx] += 1.0;
    for(i = 0; i < mot->len-3; i++) {
      for(j = start-(mot->spacer)-(mot->len); j <= start-(mot->spacer)-(i+3);
          j++) {
        if(j < 0) continue;
        if(j <= start-16-i) spacendx = 3;
        else if(j <= start-14-i) spacendx = 2;
        else if(j >= start-7-i) spacendx = 1;
        else spacendx = 0;
        mcnt[i][spacendx][mer_ndx(i+3, wseq, j)] += 1.0;
      }
    }
  }
  /* Stage 2:  Only count the highest scoring motif. */
  else if(stage == 2) mcnt[mot->len-3][mot->spacendx][mot->ndx] += 1.0;
}

/*******************************************************************************
  In addition to log likelihood, we also require a motif to actually be
  present a good portion of the time in an absolute sense across the genome.
  The coverage map is just a numerical map of whether or not to accept a
  putative motif as a real one (despite its log likelihood score.  A motif is
  considered "good" if it contains a 3-base subset of itself that is present
  in at least 20% of the total genes.  In the final stage of iterative start
  training, all motifs are labeled good.  0 = bad, 1 = good, 2 = good
  w/mismatch.
*******************************************************************************/
void build_coverage_map(double real[4][4][4096], int good[4][4][4096], double
                        ng, int stage) {
  int i, j, k, l, tmp, decomp[3];
  double thresh = 0.2;

  for(i = 0; i < 4; i++) for(j = 0; j < 4; j++) for(k = 0; k < 4096; k++)
    good[i][j][k] = 0;

  /* 3-base motifs */
  for(i = 0; i < 4; i++) for(j = 0; j < 64; j++) {
    if(real[0][i][j]/ng >= thresh) { for(k = 0; k < 4; k++) good[0][k][j] = 1; }
  }

  /* 4-base motifs, must contain two valid 3-base motifs */
  for(i = 0; i < 4; i++) for(j = 0; j < 256; j++) {
    decomp[0] = (j&252)>>2;
    decomp[1] = j&63;
    if(good[0][i][decomp[0]] == 0 || good[0][i][decomp[1]] == 0) continue;
    good[1][i][j] = 1;
  }

  /* 5-base motifs, interior mismatch allowed only if entire 5-base */
  /* motif represents 3 valid 3-base motifs (if mismatch converted) */
  for(i = 0; i < 4; i++) for(j = 0; j < 1024; j++) {
    decomp[0] = (j&1008)>>4;
    decomp[1] = (j&252)>>2;
    decomp[2] = j&63;
    if(good[0][i][decomp[0]] == 0 || good[0][i][decomp[1]] == 0 ||
       good[0][i][decomp[2]] == 0)
      continue;
    good[2][i][j] = 1;
    tmp = j;
    for(k = 0; k <= 16; k+= 16) {
      tmp = tmp ^ k;
      for(l = 0; l <= 32; l+= 32) {
        tmp = tmp ^ l;
        if(good[2][i][tmp] == 0) good[2][i][tmp] = 2;
      }
    }
  }

  /* 6-base motifs, must contain two valid 5-base motifs */
  for(i = 0; i < 4; i++) for(j = 0; j < 4096; j++) {
    decomp[0] = (j&4092)>>2;
    decomp[1] = j&1023;
    if(good[2][i][decomp[0]] == 0 || good[2][i][decomp[1]] == 0) continue;
    if(good[2][i][decomp[0]] == 1 && good[2][i][decomp[1]] == 1)
      good[3][i][j] = 1;
    else good[3][i][j] = 2;
  }

/* output all good motifs, useful info, keeping it in
printf("GOOD MOTIFS\n");
for(i = 0; i < 4; i++) for(j = 0; j < 4; j++) for(k = 0; k < 4096; k++) {
  if(good[i][j][k] == 0) continue;
  mer_text(qt, i+3, k);
  printf("motif %s %d %d %d is good\n", qt, i+3, j, k);
}
*/

}


/*******************************************************************************
  When connecting two genes, we add a bonus for the -1 and -4 base overlaps on
  the same strand, which often signify an operon and negate the need for an
  RBS for the second gene.  In addition, we add a slight bonus when genes are
  close and a slight penalty when switching strands or having a large
  intergenic space.
*******************************************************************************/
double intergenic_mod(struct _node *n1, struct _node *n2, struct _training
                      *tinf) {
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
    rval -= 0.15 * tinf->st_wt;
  }
  else if((dist <= OPER_DIST && ovlp == 0) || dist < 0.25*OPER_DIST) {
    rval += (2.0 - (double)(dist)/OPER_DIST) * 0.15 * tinf->st_wt;
  }
  return rval;
}

/*******************************************************************************
  Write detailed scoring information about every single possible gene.  Only
  done at the user's request.
*******************************************************************************/
void write_start_file(FILE *fh, struct _node *nod, int nn, struct _training
                      *tinf, int sctr, int slen, int is_meta, char *mdesc, 
                      char *version, char *header) {
  int i, prev_stop = -1, prev_strand = 0, st_type;
  double rbs1, rbs2;
  char sd_string[28][100], sd_spacer[28][20], qt[10];
  char type_string[4][5] = { "ATG", "GTG", "TTG" , "Edge" };
  char seq_data[MAX_LINE*2], run_data[MAX_LINE];
  char buffer[MAX_LINE] = {0};

  /* Initialize sequence data */
  sprintf(seq_data, "seqnum=%d;seqlen=%d;seqhdr=\"%s\"", sctr, slen, header);

  /* Initialize run data string */
  if(is_meta == 0) {
    sprintf(run_data, "version=Prodigal.v%s;run_type=Single;", version);
    strcat(run_data, "model=\"Ab initio\";");
  }
  else {
    sprintf(run_data, "version=Prodigal.v%s;run_type=Metagenomic;", version);
    sprintf(buffer, "model=\"%s\";", mdesc);
    strcat(run_data, buffer);
  }
  sprintf(buffer, "gc_cont=%.2f;transl_table=%d;uses_sd=%d",
          tinf->gc*100.0, tinf->trans_table, tinf->uses_sd);
  strcat(run_data, buffer);
 
  strcpy(sd_string[0], "None");
  strcpy(sd_spacer[0], "None");
  strcpy(sd_string[1], "GGA/GAG/AGG");
  strcpy(sd_spacer[1], "3-4bp");
  strcpy(sd_string[2], "3Base/5BMM");
  strcpy(sd_spacer[2], "13-15bp");
  strcpy(sd_string[3], "4Base/6BMM");
  strcpy(sd_spacer[3], "13-15bp");
  strcpy(sd_string[4], "AGxAG");
  strcpy(sd_spacer[4], "11-12bp");
  strcpy(sd_string[5], "AGxAG");
  strcpy(sd_spacer[5], "3-4bp");
  strcpy(sd_string[6], "GGA/GAG/AGG");
  strcpy(sd_spacer[6], "11-12bp");
  strcpy(sd_string[7], "GGxGG");
  strcpy(sd_spacer[7], "11-12bp");
  strcpy(sd_string[8], "GGxGG");
  strcpy(sd_spacer[8], "3-4bp");
  strcpy(sd_string[9], "AGxAG");
  strcpy(sd_spacer[9], "5-10bp");
  strcpy(sd_string[10], "AGGAG(G)/GGAGG");
  strcpy(sd_spacer[10], "13-15bp");
  strcpy(sd_string[11], "AGGA/GGAG/GAGG");
  strcpy(sd_spacer[11], "3-4bp");
  strcpy(sd_string[12], "AGGA/GGAG/GAGG");
  strcpy(sd_spacer[12], "11-12bp");
  strcpy(sd_string[13], "GGA/GAG/AGG");
  strcpy(sd_spacer[13], "5-10bp");
  strcpy(sd_string[14], "GGxGG");
  strcpy(sd_spacer[14], "5-10bp");
  strcpy(sd_string[15], "AGGA");
  strcpy(sd_spacer[15], "5-10bp");
  strcpy(sd_string[16], "GGAG/GAGG");
  strcpy(sd_spacer[16], "5-10bp");
  strcpy(sd_string[17], "AGxAGG/AGGxGG");
  strcpy(sd_spacer[17], "11-12bp");
  strcpy(sd_string[18], "AGxAGG/AGGxGG");
  strcpy(sd_spacer[18], "3-4bp");
  strcpy(sd_string[19], "AGxAGG/AGGxGG");
  strcpy(sd_spacer[19], "5-10bp");
  strcpy(sd_string[20], "AGGAG/GGAGG");
  strcpy(sd_spacer[20], "11-12bp");
  strcpy(sd_string[21], "AGGAG");
  strcpy(sd_spacer[21], "3-4bp");
  strcpy(sd_string[22], "AGGAG");
  strcpy(sd_spacer[22], "5-10bp");
  strcpy(sd_string[23], "GGAGG");
  strcpy(sd_spacer[23], "3-4bp");
  strcpy(sd_string[24], "GGAGG");
  strcpy(sd_spacer[24], "5-10bp");
  strcpy(sd_string[25], "AGGAGG");
  strcpy(sd_spacer[25], "11-12bp");
  strcpy(sd_string[26], "AGGAGG");
  strcpy(sd_spacer[26], "3-4bp");
  strcpy(sd_string[27], "AGGAGG");
  strcpy(sd_spacer[27], "5-10bp");

  qsort(nod, nn, sizeof(struct _node), &stopcmp_nodes);

  fprintf(fh, "# Sequence Data: %s\n", seq_data);
  fprintf(fh, "# Run Data: %s\n\n", run_data);

  fprintf(fh, "Beg\tEnd\tStd\tTotal\tCodPot\tStrtSc\tCodon\tRBSMot\t");
  fprintf(fh, "Spacer\tRBSScr\tUpsScr\tTypeScr\tGCCont\n");
  for(i = 0; i < nn; i++) {
    if(nod[i].type == STOP) continue;
    if(nod[i].edge == 1) st_type = 3;
    else st_type = nod[i].type;
    if(nod[i].stop_val != prev_stop || nod[i].strand != prev_strand) {
      prev_stop = nod[i].stop_val;
      prev_strand = nod[i].strand;
      fprintf(fh, "\n");
    }
    if(nod[i].strand == 1)
      fprintf(fh, "%d\t%d\t+\t%.2f\t%.2f\t%.2f\t%s\t", nod[i].ndx+1,
              nod[i].stop_val+3, nod[i].cscore+nod[i].sscore, nod[i].cscore,
              nod[i].sscore, type_string[st_type]);
    if(nod[i].strand == -1)
      fprintf(fh, "%d\t%d\t-\t%.2f\t%.2f\t%.2f\t%s\t", nod[i].stop_val-1,
              nod[i].ndx+1, nod[i].cscore+nod[i].sscore, nod[i].cscore,
              nod[i].sscore, type_string[st_type]);
    rbs1 = tinf->rbs_wt[nod[i].rbs[0]]*tinf->st_wt;
    rbs2 = tinf->rbs_wt[nod[i].rbs[1]]*tinf->st_wt;
    if(tinf->uses_sd == 1) {
      if(rbs1 > rbs2) {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nod[i].rbs[0]], 
                sd_spacer[nod[i].rbs[0]], nod[i].rscore);
      }
      else {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nod[i].rbs[1]], 
                sd_spacer[nod[i].rbs[1]], nod[i].rscore);
      }
    }
    else {
      mer_text(qt, nod[i].mot.len, nod[i].mot.ndx);
      if(tinf->no_mot > -0.5 && rbs1 > rbs2 && rbs1 > nod[i].mot.score *
         tinf->st_wt) {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nod[i].rbs[0]], 
                sd_spacer[nod[i].rbs[0]], nod[i].rscore);
      }
      else if(tinf->no_mot > -0.5 && rbs2 >= rbs1 && rbs2 > nod[i].mot.score *
              tinf->st_wt) {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nod[i].rbs[1]], 
                sd_spacer[nod[i].rbs[1]], nod[i].rscore);
      }
      else {
        if(nod[i].mot.len == 0) fprintf(fh, "None\tNone\t%.2f\t", 
                                      nod[i].rscore);
        else fprintf(fh, "%s\t%dbp\t%.2f\t", qt, nod[i].mot.spacer, 
                     nod[i].rscore);
      }
    }
    fprintf(fh, "%.2f\t%.2f\t%.3f\n", nod[i].uscore, nod[i].tscore, nod[i].gc_cont);
  }
  fprintf(fh, "\n");
  qsort(nod, nn, sizeof(struct _node), &compare_nodes);
}

/* Checks to see if a node boundary crosses a mask */
int cross_mask(int x, int y, mask *mlist, int nm) {
  int i;
  for(i = 0; i < nm; i++) {
    if(y < mlist[i].begin || x > mlist[i].end) continue;
    return 1;
  }
  return 0;
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
