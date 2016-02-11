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

#include "dprog.h"

/*******************************************************************************
  Basic dynamic programming routine for predicting genes.  The 'flag' variable
  is set to 0 for the initial dynamic programming routine based solely on GC
  frame plot (used to construct a training set.  If the flag is set to 1, the
  routine does the final dynamic programming based on
  coding, RBS scores, etc.
*******************************************************************************/

int dprog(struct _node *nod, int nn, struct _training *tinf, int flag) {
  int i, j, min, max_ndx = -1, path, nxt, tmp;
  double max_sc = -1.0;

  if(nn == 0) return -1;
  for(i = 0; i < nn; i++) {
    nod[i].score = 0;
    nod[i].traceb = -1;
    nod[i].tracef = -1;
  }
  for(i = 0; i < nn; i++) {

    /* Set up distance constraints for making connections, */
    /* but make exceptions for giant ORFS.                 */
    if(i < MAX_NODE_DIST) min = 0; else min = i-MAX_NODE_DIST;
    if(nod[i].strand == -1 && nod[i].type != STOP && nod[min].ndx >=
       nod[i].stop_val)
      while(min >= 0 && nod[i].ndx != nod[i].stop_val) min--;
    if(nod[i].strand == 1 && nod[i].type == STOP && nod[min].ndx >=
       nod[i].stop_val)
      while(min >= 0 && nod[i].ndx != nod[i].stop_val) min--;
    if(min < MAX_NODE_DIST) min = 0;
    else min = min-MAX_NODE_DIST;
    for(j = min; j < i; j++) {
      score_connection(nod, j, i, tinf, flag);
    }
  }
  for(i = nn-1; i >= 0; i--) {
    if(nod[i].strand == 1 && nod[i].type != STOP) continue;
    if(nod[i].strand == -1 && nod[i].type == STOP) continue;
    if(nod[i].score > max_sc) { max_sc = nod[i].score; max_ndx = i; }
  }

  /* First Pass: untangle the triple overlaps */
  path = max_ndx;
  while(nod[path].traceb != -1) {
    nxt = nod[path].traceb;
    if(nod[path].strand == -1 && nod[path].type == STOP &&
      nod[nxt].strand == 1 && nod[nxt].type == STOP &&
      nod[path].ov_mark != -1 && nod[path].ndx > nod[nxt].ndx) {
      tmp = nod[path].star_ptr[nod[path].ov_mark];
      for(i = tmp; nod[i].ndx != nod[tmp].stop_val; i--);
      nod[path].traceb = tmp;
      nod[tmp].traceb = i;
      nod[i].ov_mark = -1;
      nod[i].traceb = nxt;
    }
    path = nod[path].traceb;
  }

  /* Second Pass: Untangle the simple overlaps */
  path = max_ndx;
  while(nod[path].traceb != -1) {
    nxt = nod[path].traceb;
    if(nod[path].strand == -1 && nod[path].type != STOP && nod[nxt].strand == 1
       && nod[nxt].type == STOP) {
      for(i = path; nod[i].ndx != nod[path].stop_val; i--);
      nod[path].traceb = i; nod[i].traceb = nxt;
    }
    if(nod[path].strand == 1 && nod[path].type == STOP && nod[nxt].strand == 1
       && nod[nxt].type == STOP) {
      nod[path].traceb = nod[nxt].star_ptr[(nod[path].ndx)%3];
      nod[nod[path].traceb].traceb = nxt;
    }
    if(nod[path].strand == -1 && nod[path].type == STOP && nod[nxt].strand == -1
       && nod[nxt].type == STOP) {
      nod[path].traceb = nod[path].star_ptr[(nod[nxt].ndx)%3];
      nod[nod[path].traceb].traceb = nxt;
    }
    path = nod[path].traceb;
  }

  /* Mark forward pointers */
  path = max_ndx;
  while(nod[path].traceb != -1) {
    nod[nod[path].traceb].tracef = path;
    path = nod[path].traceb;
  }

  if(nod[max_ndx].traceb == -1) return -1;
  else return max_ndx;
}

/*******************************************************************************
  This routine scores the connection between two nodes, the most basic of which
  is 5'fwd->3'fwd (gene) and 3'rev->5'rev (rev gene).  If the connection ending
  at n2 is the maximal scoring model, it updates the pointers in the dynamic
  programming model.  n3 is used to handle overlaps, i.e. cases where 5->3'
  overlaps 5'->3' on the same strand.  In this case, 3' connects directly to 3',
  and n3 is used to untangle the 5' end of the second gene.
*******************************************************************************/

void score_connection(struct _node *nod, int p1, int p2, struct _training *tinf,
                      int flag) {
  struct _node *n1 = &(nod[p1]), *n2 = &(nod[p2]), *n3;
  int i, left = n1->ndx, right = n2->ndx, bnd, ovlp = 0, maxfr = -1;
  double score = 0.0, scr_mod = 0.0, maxval;

  /***********************/
  /* Invalid Connections */
  /***********************/

  /* 5'fwd->5'fwd, 5'rev->5'rev */
  if(n1->type != STOP && n2->type != STOP && n1->strand == n2->strand) return;

  /* 5'fwd->5'rev, 5'fwd->3'rev */
  else if(n1->strand == 1 && n1->type != STOP && n2->strand == -1) return;

  /* 3'rev->5'fwd, 3'rev->3'fwd) */
  else if(n1->strand == -1 && n1->type == STOP && n2->strand == 1) return;

  /* 5'rev->3'fwd */
  else if(n1->strand == -1 && n1->type != STOP && n2->strand == 1 && n2->type ==
          STOP) return;

  /******************/
  /* Edge Artifacts */
  /******************/
  if(n1->traceb == -1 && n1->strand == 1 && n1->type == STOP) return;
  if(n1->traceb == -1 && n1->strand == -1 && n1->type != STOP) return;

  /*********/
  /* Genes */
  /*********/

  /* 5'fwd->3'fwd */
  else if(n1->strand == n2->strand && n1->strand == 1 && n1->type != STOP &&
          n2->type == STOP) {
    if(n2->stop_val >= n1->ndx) return;
    if(n1->ndx % 3 != n2->ndx % 3) return;
    right += 2;
    if(flag == 0) scr_mod = tinf->bias[0]*n1->gc_score[0] +
                  tinf->bias[1]*n1->gc_score[1] + tinf->bias[2]*n1->gc_score[2];
    else if(flag == 1) score = n1->cscore + n1->sscore;
  }

  /* 3'rev->5'rev */
  else if(n1->strand == n2->strand && n1->strand == -1 && n1->type == STOP &&
          n2->type != STOP) {
    if(n1->stop_val <= n2->ndx) return;
    if(n1->ndx % 3 != n2->ndx % 3) return;
    left -= 2;
    if(flag == 0) scr_mod = tinf->bias[0]*n2->gc_score[0] +
                  tinf->bias[1]*n2->gc_score[1] + tinf->bias[2]*n2->gc_score[2];
    else if(flag == 1) score = n2->cscore + n2->sscore;
  }

  /********************************/
  /* Intergenic Space (Noncoding) */
  /********************************/

  /* 3'fwd->5'fwd */
  else if(n1->strand == 1 && n1->type == STOP && n2->strand == 1 && n2->type !=
          STOP) {
    left += 2;
    if(left >= right) return;
    if(flag == 1) score = intergenic_mod(n1, n2, tinf);
  }

  /* 3'fwd->3'rev */
  else if(n1->strand == 1 && n1->type == STOP && n2->strand == -1 && n2->type ==
          STOP) {
    left += 2; right -= 2;
    if(left >= right) return;
    /* Overlapping Gene Case 2: Three consecutive overlapping genes f r r */
    maxfr = -1; maxval = 0.0;
    for(i = 0; i < 3; i++) {
      if(n2->star_ptr[i] == -1) continue;
      n3 = &(nod[n2->star_ptr[i]]);
      ovlp = left - n3->stop_val + 3;
      if(ovlp <= 0 || ovlp >= MAX_OPP_OVLP) continue;
      if(ovlp >= n3->ndx - left) continue;
      if(n1->traceb == -1) continue;
      if(ovlp >= n3->stop_val - nod[n1->traceb].ndx - 2) continue;
      if((flag == 1 && n3->cscore + n3->sscore + intergenic_mod(n3, n2, tinf) >
          maxval) || (flag == 0 && tinf->bias[0]*n3->gc_score[0] +
          tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2] >
          maxval)) {
        maxfr = i;
        maxval = n3->cscore + n3->sscore + intergenic_mod(n3, n2, tinf);
      }
    }
    if(maxfr != -1) {
      n3 = &(nod[n2->star_ptr[maxfr]]);
      if(flag == 0) scr_mod = tinf->bias[0]*n3->gc_score[0] +
                    tinf->bias[1]*n3->gc_score[1] +
                    tinf->bias[2]*n3->gc_score[2];
      else if(flag == 1) score = n3->cscore + n3->sscore +
              intergenic_mod(n3, n2, tinf);
    }
    else if(flag == 1) score = intergenic_mod(n1, n2, tinf);
  }

  /* 5'rev->3'rev */
  else if(n1->strand == -1 && n1->type != STOP && n2->strand == -1 && n2->type
          == STOP) {
    right -= 2;
    if(left >= right) return;
    if(flag == 1) score = intergenic_mod(n1, n2, tinf);
  }

  /* 5'rev->5'fwd */
  else if(n1->strand == -1 && n1->type != STOP && n2->strand == 1 && n2->type
          != STOP) {
    if(left >= right) return;
    if(flag == 1) score = intergenic_mod(n1, n2, tinf);
  }

  /********************/
  /* Possible Operons */
  /********************/

  /* 3'fwd->3'fwd, check for a start just to left of first 3' */
  else if(n1->strand == 1 && n2->strand == 1 && n1->type == STOP && n2->type ==
          STOP) {
    if(n2->stop_val >= n1->ndx) return;
    if(n1->star_ptr[n2->ndx%3] == -1) return;
    n3 = &(nod[n1->star_ptr[n2->ndx%3]]);
    left = n3->ndx; right += 2;
    if(flag == 0) scr_mod = tinf->bias[0]*n3->gc_score[0] +
                  tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2];
    else if(flag == 1) score = n3->cscore + n3->sscore +
                       intergenic_mod(n1, n3, tinf);
  }

  /* 3'rev->3'rev, check for a start just to right of second 3' */
  else if(n1->strand == -1 && n1->type == STOP && n2->strand == -1 && n2->type
          == STOP) {
    if(n1->stop_val <= n2->ndx) return;
    if(n2->star_ptr[n1->ndx%3] == -1) return;
    n3 = &(nod[n2->star_ptr[n1->ndx%3]]);
    left -= 2; right = n3->ndx;
    if(flag == 0) scr_mod = tinf->bias[0]*n3->gc_score[0] +
                  tinf->bias[1]*n3->gc_score[1] + tinf->bias[2]*n3->gc_score[2];
    else if(flag == 1) score = n3->cscore + n3->sscore +
                       intergenic_mod(n3, n2, tinf);
  }

  /***************************************/
  /* Overlapping Opposite Strand 3' Ends */
  /***************************************/

  /* 3'for->5'rev */
  else if(n1->strand == 1 && n1->type == STOP && n2->strand == -1 && n2->type
          != STOP) {
    if(n2->stop_val-2 >= n1->ndx+2) return;
    ovlp = (n1->ndx+2) - (n2->stop_val-2) + 1;
    if(ovlp >= MAX_OPP_OVLP) return;
    if((n1->ndx+2 - n2->stop_val-2 + 1) >= (n2->ndx -n1->ndx+3 + 1)) return;
    if(n1->traceb == -1) bnd = 0;
    else bnd = nod[n1->traceb].ndx;
    if((n1->ndx+2 - n2->stop_val-2 + 1) >= (n2->stop_val-3 - bnd + 1)) return;
    left = n2->stop_val-2;
    if(flag == 0) scr_mod = tinf->bias[0]*n2->gc_score[0] +
                  tinf->bias[1]*n2->gc_score[1] + tinf->bias[2]*n2->gc_score[2];
    else if(flag == 1) score = n2->cscore + n2->sscore - 0.15*tinf->st_wt;
  }

  if(flag == 0) score = ((double)(right-left+1-(ovlp*2)))*scr_mod;

  if(n1->score + score >= n2->score) {
    n2->score = n1->score + score;
    n2->traceb = p1;
    n2->ov_mark = maxfr;
  }

  return;
}

/*******************************************************************************
  Sometimes bad genes creep into the model due to the node distance constraint
  in the dynamic programming routine.  This routine just does a sweep through
  the genes and eliminates ones with negative scores.
*******************************************************************************/

void eliminate_bad_genes(struct _node *nod, int dbeg, struct _training *tinf) {
  int path;

  if(dbeg == -1) return;
  path = dbeg;
  while(nod[path].traceb != -1) path = nod[path].traceb;
  while(nod[path].tracef != -1) {
    if(nod[path].strand == 1 && nod[path].type == STOP)
      nod[nod[path].tracef].sscore += intergenic_mod(&nod[path],
                                      &nod[nod[path].tracef], tinf);
    if(nod[path].strand == -1 && nod[path].type != STOP)
      nod[path].sscore += intergenic_mod(&nod[path], &nod[nod[path].tracef],
                          tinf);
    path = nod[path].tracef;
  }

  path = dbeg;
  while(nod[path].traceb != -1) path = nod[path].traceb;
  while(nod[path].tracef != -1) {
    if(nod[path].strand == 1 && nod[path].type != STOP &&
       nod[path].cscore + nod[path].sscore < 0) {
      nod[path].elim = 1; nod[nod[path].tracef].elim = 1;
    }
    if(nod[path].strand == -1 && nod[path].type == STOP &&
       nod[nod[path].tracef].cscore + nod[nod[path].tracef].sscore < 0) {
      nod[path].elim = 1; nod[nod[path].tracef].elim = 1;
    }
    path = nod[path].tracef;
  }
}
