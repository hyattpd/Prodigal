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

#include "gene.h"

/* Copies genes from the dynamic programming to a final array */

int add_genes(struct _gene *glist, struct _node *nod, int dbeg) {
  int path, ctr;

  if(dbeg == -1) return 0;
  path = dbeg; ctr = 0;
  while(nod[path].traceb != -1) path = nod[path].traceb;

  while(path != -1) {
    if(nod[path].elim == 1) { path = nod[path].tracef; continue; }
    if(nod[path].strand == 1 && nod[path].type != STOP) {
      glist[ctr].begin = nod[path].ndx+1;
      glist[ctr].start_ndx = path;
    }
    if(nod[path].strand == -1 && nod[path].type == STOP) {
      glist[ctr].begin = nod[path].ndx-1;
      glist[ctr].stop_ndx = path;
    }
    if(nod[path].strand == 1 && nod[path].type == STOP) {
      glist[ctr].end = nod[path].ndx+3;
      glist[ctr].stop_ndx = path;
      ctr++;
    }
    if(nod[path].strand == -1 && nod[path].type != STOP) {
      glist[ctr].end = nod[path].ndx+1;
      glist[ctr].start_ndx = path;
      ctr++;
    }
    path = nod[path].tracef;
    if(ctr == MAX_GENES) {
#ifndef PRODIGAL_NO_MAIN
      fprintf(stderr, "warning, max # of genes exceeded, truncating...\n");
#endif
      return ctr;
    }
  }
  return ctr;
}

/*******************************************************************************
  This routine attempts to solve the problem of extremely close starts.  If two
  potential starts are 5 amino acids or less away from each other, this routine
  sets their coding equal to each other and lets the RBS/operon/ATG-GTG-TTG
  start features determine which start to use, under the assumption that 5 or
  less words of coding is too weak a signal to use to select the proper start.

  In addition, we try to correct TTG (or whatever start codon is rare) starts
  that have an RBS score, an upstream score, and a coding score all superior
  to whatever start we initially chose.

  This routine was tested on numerous genomes and found to increase overall
  performance. 
*******************************************************************************/
void tweak_final_starts(struct _gene *genes, int ng, struct _node *nod,
                        int nn, struct _training *tinf) {
  int i, j, ndx, mndx, maxndx[2];
  double sc, igm, tigm, maxsc[2], maxigm[2];

  for(i = 0; i < ng; i++) {
    ndx = genes[i].start_ndx;
    sc = nod[ndx].sscore + nod[ndx].cscore;
    igm = 0.0;
    if(i > 0 && nod[ndx].strand == 1 && nod[genes[i-1].start_ndx].strand == 1)
      igm = intergenic_mod(&nod[genes[i-1].stop_ndx], &nod[ndx], tinf);
    if(i > 0 && nod[ndx].strand == 1 && nod[genes[i-1].start_ndx].strand == -1)
      igm = intergenic_mod(&nod[genes[i-1].start_ndx], &nod[ndx], tinf);
    if(i < ng-1 && nod[ndx].strand == -1 && nod[genes[i+1].start_ndx].strand 
       == 1)
      igm = intergenic_mod(&nod[ndx], &nod[genes[i+1].start_ndx], tinf);
    if(i < ng-1 && nod[ndx].strand == -1 && nod[genes[i+1].start_ndx].strand 
       == -1)
      igm = intergenic_mod(&nod[ndx], &nod[genes[i+1].stop_ndx], tinf);

    /* Search upstream and downstream for the #2 and #3 scoring starts */
    maxndx[0] = -1; maxndx[1] = -1; maxsc[0] = 0; maxsc[1] = 0;
    maxigm[0] = 0; maxigm[1] = 0;
    for(j = ndx-100; j < ndx+100; j++) {
      if(j < 0 || j >= nn || j == ndx) continue;
      if(nod[j].type == STOP || nod[j].stop_val != nod[ndx].stop_val)
        continue;

      tigm = 0.0;
      if(i > 0 && nod[j].strand == 1 && nod[genes[i-1].start_ndx].strand == 1)
      {
        if(nod[genes[i-1].stop_ndx].ndx - nod[j].ndx > MAX_SAM_OVLP) continue;
        tigm = intergenic_mod(&nod[genes[i-1].stop_ndx], &nod[j], tinf);
      }
      if(i > 0 && nod[j].strand == 1 && nod[genes[i-1].start_ndx].strand == -1)
      {
        if(nod[genes[i-1].start_ndx].ndx - nod[j].ndx >= 0) continue;
        tigm = intergenic_mod(&nod[genes[i-1].start_ndx], &nod[j], tinf);
      }
      if(i < ng-1 && nod[j].strand == -1 && nod[genes[i+1].start_ndx].strand 
         == 1) {
        if(nod[j].ndx - nod[genes[i+1].start_ndx].ndx >= 0) continue;
        tigm = intergenic_mod(&nod[j], &nod[genes[i+1].start_ndx], tinf);
      }
      if(i < ng-1 && nod[j].strand == -1 && nod[genes[i+1].start_ndx].strand 
         == -1) {
        if(nod[j].ndx - nod[genes[i+1].stop_ndx].ndx > MAX_SAM_OVLP) continue;
        tigm = intergenic_mod(&nod[j], &nod[genes[i+1].stop_ndx], tinf);
      }
 
      if(maxndx[0] == -1) {
        maxndx[0] = j;
        maxsc[0] = nod[j].cscore + nod[j].sscore;
        maxigm[0] = tigm;
      }
      else if(nod[j].cscore + nod[j].sscore + tigm > maxsc[0]) {
        maxndx[1] = maxndx[0];
        maxsc[1] = maxsc[0];
        maxigm[1] = maxigm[0];
        maxndx[0] = j;
        maxsc[0] = nod[j].cscore + nod[j].sscore;
        maxigm[0] = tigm;
      }
      else if(maxndx[1] == -1 || nod[j].cscore + nod[j].sscore + tigm > 
              maxsc[1]) { 
        maxndx[1] = j;
        maxsc[1] = nod[j].cscore + nod[j].sscore;
        maxigm[1] = tigm;
      }
    }

    /* Change the start if it's a TTG with better coding/RBS/upstream score */
    /* Also change the start if it's <=15bp but has better coding/RBS       */
    for(j = 0; j < 2; j++) {
      mndx = maxndx[j];
      if(mndx == -1) continue;

      /* Start of less common type but with better coding, rbs, and */
      /* upstream.  Must be 18 or more bases away from original.    */
      if(nod[mndx].tscore < nod[ndx].tscore && maxsc[j]-nod[mndx].tscore >= 
         sc-nod[ndx].tscore+tinf->st_wt && nod[mndx].rscore > nod[ndx].rscore
         && nod[mndx].uscore > nod[ndx].uscore && nod[mndx].cscore > 
         nod[ndx].cscore && abs(nod[mndx].ndx-nod[ndx].ndx) > 15) {
        maxsc[j] += nod[ndx].tscore-nod[mndx].tscore;
      }

      /* Close starts.  Ignore coding and see if start has better rbs */
      /* and type. */
      else if(abs(nod[mndx].ndx-nod[ndx].ndx) <= 15 && nod[mndx].rscore+
              nod[mndx].tscore > nod[ndx].rscore+nod[ndx].tscore &&
              nod[ndx].edge == 0 && nod[mndx].edge == 0) {
        if(nod[ndx].cscore > nod[mndx].cscore) 
          maxsc[j] += nod[ndx].cscore - nod[mndx].cscore;
        if(nod[ndx].uscore > nod[mndx].uscore) 
          maxsc[j] += nod[ndx].uscore - nod[mndx].uscore;
        if(igm > maxigm[j]) maxsc[j] += igm - maxigm[j]; 
      }
  
      else maxsc[j] = -1000.0;
    }

    /* Change the gene coordinates to the new maximum. */
    mndx = -1;
    for(j = 0; j < 2; j++) {
      if(maxndx[j] == -1) continue;
      if(mndx == -1 && maxsc[j]+maxigm[j] > sc+igm)
        mndx = j;
      else if(mndx >= 0 && maxsc[j]+maxigm[j] > maxsc[mndx]+maxigm[mndx])
        mndx = j; 
    }
    if(mndx != -1 && nod[maxndx[mndx]].strand == 1) {
      genes[i].start_ndx = maxndx[mndx];
      genes[i].begin = nod[maxndx[mndx]].ndx+1;
    } 
    else if(mndx != -1 && nod[maxndx[mndx]].strand == -1) {
      genes[i].start_ndx = maxndx[mndx];
      genes[i].end = nod[maxndx[mndx]].ndx+1;
    } 
  }
}

void record_gene_data(struct _gene *genes, int ng, struct _node *nod,
                      struct _training *tinf, int sctr) {
  int i;
  char type_string[4][5] = { "ATG", "GTG", "TTG" , "Edge" };

  for(i = 0; i < ng; i++) {
    struct _gene *g = &genes[i];
    struct _node *st = &nod[g->start_ndx], *sp = &nod[g->stop_ndx];

    /* Record basic gene data */
    g->strand = st->strand;
    g->partial_left = (st->edge == 1 && st->strand == 1) ||
                      (sp->edge == 1 && st->strand == -1);
    g->partial_right = (sp->edge == 1 && st->strand == 1) ||
                       (st->edge == 1 && st->strand == -1);
    g->start_type = (st->edge == 1) ? 3 : st->type;
    pdg_rbs_text(&g->rbs, st, tinf);
    g->gc_cont = st->gc_cont;

    /* Record score data */
    g->cscore = st->cscore;
    g->sscore = st->sscore;
    g->rscore = st->rscore;
    g->uscore = st->uscore;
    g->tscore = st->tscore;
    g->confidence = calculate_confidence(st->cscore + st->sscore, tinf->st_wt);

    sprintf(g->gene_data, "ID=%d_%d;partial=%d%d;start_type=%s;"
            "rbs_motif=%s;rbs_spacer=%s;gc_cont=%.3f", sctr, i+1,
            g->partial_left, g->partial_right, type_string[g->start_type],
            g->rbs.motif, g->rbs.spacer, g->gc_cont);
    sprintf(g->score_data, "conf=%.2f;score=%.2f;cscore=%.2f;sscore=%.2f;"
            "rscore=%.2f;uscore=%.2f;tscore=%.2f;", g->confidence,
            g->cscore + g->sscore, g->cscore, g->sscore, g->rscore,
            g->uscore, g->tscore);
  }
}

/* Print the genes.  'Flag' indicates which format to use. */
void print_genes(FILE *fp, struct _gene *genes, int ng, struct _node *nod, 
                 int slen, int flag, int sctr, int is_meta, char *mdesc,
                 struct _training *tinf, char *header, char *short_hdr,
                 char *version) {
  int i, ndx, sndx;
  char left[50], right[50];
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

  strcpy(left, "");
  strcpy(right, "");

  /* Print the gff header once */
  if(flag == 3 && sctr == 1) fprintf(fp, "##gff-version  3\n");

  /* Print sequence/model information */
  if(flag == 0) {
    fprintf(fp, "DEFINITION  %s;%s\n", seq_data, run_data);
    fprintf(fp, "FEATURES             Location/Qualifiers\n");
  }
  else if(flag != 1) {
    fprintf(fp, "# Sequence Data: %s\n", seq_data);
    fprintf(fp, "# Model Data: %s\n", run_data);
  }
  
  /* Print the genes */
  for(i = 0; i < ng; i++) {
    ndx = genes[i].start_ndx;
    sndx = genes[i].stop_ndx;

    /* Print the coordinates and data */
    if(nod[ndx].strand == 1) {

      if(nod[ndx].edge == 1) sprintf(left, "<%d", genes[i].begin);
      else sprintf(left, "%d", genes[i].begin);
      if(nod[sndx].edge == 1) sprintf(right, ">%d", genes[i].end);
      else sprintf(right, "%d", genes[i].end);

      if(flag == 0) {
        fprintf(fp, "     CDS             %s..%s\n", left, right);
        fprintf(fp, "                     ");
        fprintf(fp, "/note=\"%s;%s\"\n", genes[i].gene_data,
                genes[i].score_data);
      }
      if(flag == 1)
        fprintf(fp, "gene_prodigal=%d|1|f|y|y|3|0|%d|%d|%d|%d|-1|-1|1.0\n", i+1,
                genes[i].begin, genes[i].end, genes[i].begin, genes[i].end);
      if(flag == 2) fprintf(fp, ">%d_%d_%d_+\n", i+1, genes[i].begin, 
                            genes[i].end);
      if(flag == 3) {
        fprintf(fp, "%s\tProdigal_v%s\tCDS\t%d\t%d\t%.1f\t+\t0\t%s;%s", 
                short_hdr, version, genes[i].begin, genes[i].end, 
                nod[ndx].cscore+nod[ndx].sscore, genes[i].gene_data,
                genes[i].score_data);
        fprintf(fp, "\n"); 
      }
    }
    else {

      if(nod[sndx].edge == 1) sprintf(left, "<%d", genes[i].begin);
      else sprintf(left, "%d", genes[i].begin);
      if(nod[ndx].edge == 1) sprintf(right, ">%d", genes[i].end);
      else sprintf(right, "%d", genes[i].end);

      if(flag == 0) {
        fprintf(fp, "     CDS             complement(%s..%s)\n", left, right);
        fprintf(fp, "                     ");
        fprintf(fp, "/note=\"%s;%s\"\n", genes[i].gene_data,
                genes[i].score_data);
      }
      if(flag == 1)
        fprintf(fp, "gene_prodigal=%d|1|r|y|y|3|0|%d|%d|%d|%d|-1|-1|1.0\n", i+1,
               slen+1-genes[i].end, slen+1-genes[i].begin,
               slen+1-genes[i].end, slen+1-genes[i].begin);
      if(flag == 2) fprintf(fp, ">%d_%d_%d_-\n", i+1, genes[i].begin, 
                            genes[i].end);
      if(flag == 3) {
        fprintf(fp, "%s\tProdigal_v%s\tCDS\t%d\t%d\t%.1f\t-\t0\t%s;%s",
                short_hdr, version, genes[i].begin, genes[i].end, 
                nod[ndx].cscore+nod[ndx].sscore, genes[i].gene_data,
                genes[i].score_data);
        fprintf(fp, "\n"); 
      }
    }
  }

  /* Footer */
  if(flag == 0) fprintf(fp, "//\n");

}

/* Print the gene translations */
void write_translations(FILE *fh, struct _gene *genes, int ng, struct 
                        _node *nod, unsigned char *seq, unsigned char *rseq, 
                        unsigned char *useq, int slen, struct _training *tinf,
                        int sctr, char *short_hdr) {
  int i, j;

  for(i = 0; i < ng; i++) {
    if(nod[genes[i].start_ndx].strand == 1) {
      fprintf(fh, ">%s_%d # %d # %d # 1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, genes[i].gene_data);
      for(j = genes[i].begin; j < genes[i].end; j+=3) {
        if(is_n(useq, j-1) == 1 || is_n(useq, j) == 1 || is_n(useq, j+1) == 1) 
          fprintf(fh, "X");
        else fprintf(fh, "%c", amino(seq, j-1, tinf, (j==genes[i].begin?1:0) &&
                     (1-nod[genes[i].start_ndx].edge)));
        if((j-genes[i].begin)%180 == 177) fprintf(fh, "\n");
      }
      if((j-genes[i].begin)%180 != 0) fprintf(fh, "\n");
    }
    else {
      fprintf(fh, ">%s_%d # %d # %d # -1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, genes[i].gene_data);
      for(j = slen+1-genes[i].end; j < slen+1-genes[i].begin; j+=3) {
        if(is_n(useq, slen-j) == 1 || is_n(useq, slen-1-j) == 1 ||
           is_n(useq, slen-2-j) == 1)
          fprintf(fh, "X");
        else fprintf(fh, "%c", amino(rseq, j-1, tinf, (j==slen+1-genes[i].end?1:0)
                     && (1-nod[genes[i].start_ndx].edge)));
        if((j-slen-1+genes[i].end)%180 == 177) fprintf(fh, "\n");
      }
      if((j-slen-1+genes[i].end)%180 != 0) fprintf(fh, "\n");
    }
  }
}

/* Print the gene nucleotide sequences */
void write_nucleotide_seqs(FILE *fh, struct _gene *genes, int ng, struct 
                           _node *nod, unsigned char *seq, unsigned char *rseq,
                           unsigned char *useq, int slen, struct _training 
                           *tinf, int sctr, char *short_hdr) {
  int i, j;

  for(i = 0; i < ng; i++) {
    if(nod[genes[i].start_ndx].strand == 1) {
      fprintf(fh, ">%s_%d # %d # %d # 1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, genes[i].gene_data);
      for(j = genes[i].begin-1; j < genes[i].end; j++) {
        if(is_a(seq, j) == 1) fprintf(fh, "A");
        else if(is_t(seq, j) == 1) fprintf(fh, "T");
        else if(is_g(seq, j) == 1) fprintf(fh, "G");
        else if(is_c(seq, j) == 1 && is_n(useq, j) == 0) fprintf(fh, "C");
        else fprintf(fh, "N");
        if((j-genes[i].begin+1)%70 == 69) fprintf(fh, "\n");
      }
      if((j-genes[i].begin+1)%70 != 0) fprintf(fh, "\n");
    }
    else {
      fprintf(fh, ">%s_%d # %d # %d # -1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, genes[i].gene_data);
      for(j = slen-genes[i].end; j < slen+1-genes[i].begin; j++) {
        if(is_a(rseq, j) == 1) fprintf(fh, "A");
        else if(is_t(rseq, j) == 1) fprintf(fh, "T");
        else if(is_g(rseq, j) == 1) fprintf(fh, "G");
        else if(is_c(rseq, j) == 1 && is_n(useq, slen-1-j) == 0) 
          fprintf(fh, "C");
        else fprintf(fh, "N");
        if((j-slen+genes[i].end)%70 == 69) fprintf(fh, "\n");
      }
      if((j-slen+genes[i].end)%70 != 0) fprintf(fh, "\n");
    }
  }
}

/* Convert score to a percent confidence */
double calculate_confidence(double score, double start_weight) {
  double conf;

  if(score/start_weight < 41) {
    conf = exp(score/start_weight);
    conf = (conf/(conf+1))*100.0;
  }
  else conf = 99.99;
  if(conf <= 50.00) { conf = 50.00; }
  return conf;
}
