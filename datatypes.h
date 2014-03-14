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

#ifndef _DATA_TYPES_H
#define _DATA_TYPES_H

struct _motif {
  int ndx;             /* Index of the best motif for this node */
  int len;             /* Length of the motif */
  int spacer;          /* Spacer between coding start and the motif */
  int spacendx;        /* Index for this spacer length */
  double score;        /* Score for the motif */
};

struct _node {
  int type;            /* 0=ATG, 1=GTG, 2=TTG/Other, 3=Stop */
  int edge;            /* Runs off the edge; 0 = normal, 1 = edge node */
  int ndx;             /* position in the sequence of the node */
  int strand;          /* 1 = forward, -1 = reverse */
  int stop_val;        /* For a stop, record previous stop; for start, record
                          its stop */
  int star_ptr[3];     /* Array of starts w/in MAX_SAM_OVLP bases of a stop in 3
                          frames */
  int gc_bias;         /* Frame of highest GC content within this node */
  double gc_score[3];  /* % GC content in different codon positions */
  double cscore;       /* Coding score for this node (based on 6-mer usage) */
  double gc_cont;      /* GC Content for the node */
  int rbs[2];          /* SD RBS score for this node (based on binding energy)
                          rbs[0] = best motif with exact match, rbs[1] = with
                          mismatches */
  struct _motif mot;   /* Upstream motif information for this node */
  double uscore;       /* Score for the upstream -1/-2, -15to-45 region */
  double tscore;       /* Score for the ATG/GTG/TTG value */
  double rscore;       /* Score for the RBS motif */
  double sscore;       /* Score for the strength of the start codon */
  int traceb;          /* Traceback to connecting node */
  int tracef;          /* Forward trace */
  int ov_mark;         /* Marker to help untangle overlapping genes */
  double score;        /* Score of total solution to this point */
  int elim;            /* If set to 1, eliminate this gene from the model */
};

struct _training {
  double gc;                    /* GC Content */
  int trans_table;              /* 11 = Standard Microbial, NCBI Trans Table to
                                   use */
  double st_wt;                 /* Start weight */
  double bias[3];               /* GC frame bias for each of the 3 positions */
  double type_wt[3];            /* Weights for ATG vs GTG vs TTG */
  int uses_sd;                  /* 0 if doesn't use SD motif, 1 if it does */
  double rbs_wt[28];            /* Set of weights for RBS scores */
  double ups_comp[32][4];       /* Base composition weights for non-RBS-distance
                                   motifs.  0-1 are the -1/-2 position, 2-31 are
                                   the -15 to -44 positions.  Second array is
                                   the base A,C,T,G,etc. */
  double mot_wt[4][4][4096];    /* Weights for upstream motifs.  First index is
                                   the motif length (3-6), the second is the
                                   spacer distance (0 = 5-10bp, 1 = 3-4bp, 2 =
                                   11-12bp, 3 = 13-15bp), and the last is the
                                   numerical value of the motif (ranging from 0
                                   to 4095 for 6-mers, less for shorter
                                   motifs) */
  double no_mot;                /* Weight for the case of no motif */
  double gene_dc[4096];         /* Coding statistics for the genome */
};

#endif
