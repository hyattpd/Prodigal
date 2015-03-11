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

#ifndef _DATA_TYPES_H
#define _DATA_TYPES_H

/* Basic structure for motifs */
struct _motif
{
  int index;            /* Index of the best motif for this node */
  int len;              /* Length of the motif */
  int spacer;           /* Spacer between coding start and the motif */
  int spacer_index;     /* Index for this spacer length */
  double score;         /* Score for the motif */
};

/* Dynamic programming node, where each node is either a start or a stop */
struct _node
{
  int type;            /* 0 = Start node, 1 = Stop node */
  int subtype;         /* For start node, 0=ATG,1=GTG,2=TTG,3=Nonst.,4=Edge */
                       /* For stop node, 0=TAA,1=TAG,2=TGA,3=Edge */
  int edge;            /* Runs off the edge: 0 = normal, 1 = edge node */
  int index;           /* position in the sequence of the node */
  int strand;          /* 1 = forward, -1 = reverse */
  int stop_val;        /* For a stop, record previous stop; for a start, record
                          its stop */
  int stop_type;       /* For starts: record subtype of matching stop node */
  int start_ptr[3];    /* Array of starts w/in MAX_SAM_OVLP bases of stop in 3
                          frames */
  int context[45];     /* Base composition around start */
  double cscore;       /* Coding score for this node */
  double lscore;       /* Length factor to the score */
  double gc_cont;      /* GC Content for the node */
  int gc_frame[3];     /* Counter for best GC frames throughout the node */
  int rbs[2];          /* SD RBS score for this node (based on binding energy)
                          rbs[0] = best motif with exact match, rbs[1] = with
                          mismatches */
  struct _motif mot;   /* Upstream motif information for this node */
  double bscore;       /* Score for the sequence context around start */
  double tscore;       /* Score for the ATG/GTG/TTG value */
  double rscore;       /* Score for the RBS motif */
  double sscore;       /* Score for the strength of the start codon */
  int trace_back;      /* Traceback to connecting node */
  int trace_forward;   /* Forward trace */
  int overlap_frame;   /* Points to best frame containing overlapping gene */
  double score;        /* Score of total solution to this point */
  int status;          /* 0 = not in the model, 1 = in the predictions */
};

/* Attributes we've studied/learned for a particular genome */
struct _training
{
  double gc;                   /* GC Content */
  int trans_table;             /* 11 = Standard Microbial, NCBI Trans Table to
                                  use */
  int uses_sd;                 /* 0 = doesn't use SD motif, 1 = it does */
  double start_weight;         /* Start weight */
  double prob_start;           /* Base probability a codon is a start */
  double prob_stop;            /* Base probability a codon is a stop */
  double type_wt[4];           /* Weights for ATG vs GTG vs TTG */
  double stop_wt[4];           /* Weights for TAA vs TGA vs TAG */
  double rbs_wt[28];           /* Set of weights for RBS scores */
  double context_wt[45][4];    /* Base composition weights for non-RBS-distance
                                  motifs.  0-1 are the -1/-2 position, 2-31 are
                                  the -15 to -44 positions.  Second array is
                                  the base A,C,T,G,etc. */
  double mot_wt[4][4][4096];   /* Weights for upstream motifs.  First index is
                                  the motif length (3-6), the second is the
                                  spacer distance (0 = 5-10bp, 1 = 3-4bp, 2 =
                                  11-12bp, 3 = 13-15bp), and the last is the
                                  numerical value of the motif (ranging from 0
                                  to 4095 for 6-mers, less with shorter
                                  motifs) */
  double no_mot;               /* Weight for the case of no motif */
  double gene_dc[4096];        /* Coding statistics for the genome */
};

/* Struct for genes, with start/stop and link to node structure */
struct _gene
{
  int begin;               /* Left end of the gene */
  int end;                 /* Right end of the gene */
  int start_index;         /* Index to the start node in the nodes file */
  int stop_index;          /* Index to the stop node in the nodes file */
  double intergenic_score; /* Score for the intergenic distances */
};

/* Structure to store text about a particular gene */
struct _gene_data
{
  char gene_data[500];     /* String containing gene information */
  char score_data[500];    /* String containing scoring information */
};

/* Summary statistics for an entire run */
struct _summary
{
  int num_contig;
  double avg_contig_len;
  int contig_len_bins[10]; /* log base 10 scale */
  double avg_contig_gc;
  int contig_gc_bins[12]; /* <25, 25-30, etc., >=75 */
  int anon_bins[50];

  int num_complete_genes;
  double avg_comp_gene_len;
  int comp_gene_len_bins[26]; /* <200, 200-400, >=5000 */
  double avg_comp_gene_gc;
  int comp_gene_gc_bins[12];
  int comp_start_bins[5]; /* ATG, GTG, TTG, Edge, Nonstandard */
  int comp_stop_bins[4]; /* TAA, TAG, TGA, Edge */
  int comp_rbs_nosd_bins[4][4][4096];

  int num_partial_genes;
  double avg_part_gene_len;
  int part_gene_len_bins[10];
  double avg_part_gene_gc;
  int part_gene_gc_bins[12];
  int part_start_bins[5];
  int part_stop_bins[4];
  int part_rbs_nosd_bins[4][4][4096];
};

/* Preset training file (with metadata) used for anonymous runs */
struct _preset_genome_bin
{
  int index;                    /* Index used for sorting */
  int clusnum;                  /* Cluster number */
  char desc[500];               /* Text description of this bin */
  double weight;                /* Current weight/score of this bin */
  double gc;                    /* GC distance from target sequence */
  struct _training *data;       /* Pointer to the training file for this bin */
};

#endif
