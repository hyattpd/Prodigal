/******************************************************************************
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
******************************************************************************/

#include "gene.h"

/******************************************************************************
  Copies genes from the dynamic programming to a more normal structure
  with starts and stops labeled, etc.  Returns the number of genes we found.
******************************************************************************/
int add_genes(struct _gene *genes, struct _node *nodes, int initial_node)
{
  int path = initial_node;   /* Traversal variable */
  int counter = 0;           /* Count number of genes */

  /* If there are no genes, return */
  if (initial_node == -1)
  {
    return 0;
  }
  /* Find the genes */
  while (path != -1)
  {
    /* Ignore deleted genes */
    if (nodes[path].status == 0)
    {
      path = nodes[path].trace_forward;
      continue;
    }
    /* Valid start/stop: add it to the list */
    if (nodes[path].strand == 1 && nodes[path].type == START)
    {
      genes[counter].begin = nodes[path].index+1;
      genes[counter].start_index = path;
    }
    else if (nodes[path].strand == -1 && nodes[path].type == STOP)
    {
      genes[counter].begin = nodes[path].index-1;
      genes[counter].stop_index = path;
    }
    else if (nodes[path].strand == 1 && nodes[path].type == STOP)
    {
      genes[counter].end = nodes[path].index+3;
      genes[counter].stop_index = path;
      counter++;
    }
    else if (nodes[path].strand == -1 && nodes[path].type == START)
    {
      genes[counter].end = nodes[path].index+1;
      genes[counter].start_index = path;
      counter++;
    }
    path = nodes[path].trace_forward;
    /* If we see too many genes (unlikely), stop adding them */
    if (counter == MAX_GENES)
    {
      fprintf(stderr, "warning, max # of genes exceeded, truncating...\n");
      return counter;
    }
  }
  return counter;
}

/******************************************************************************
  This routine attempts to solve the problem of extremely close starts.  If two
  potential starts are 5 amino acids or less away from each other, this routine
  sets their coding equal to each other and lets the RBS/operon/ATG-GTG-TTG
  start features determine which start to use, under the assumption that 5 or
  less words of coding is too weak a signal to use to select the proper start.

  In addition, we try to correct TTG (or whatever start codon is rare) starts
  that have an RBS score, an upstream score, and a coding score all superior
  to whatever start we initially chose.

  This routine was tested on numerous genomes and found to increase overall
  performance.  Having said that, this function is kind of clunky and we
  should really be able to fix the dynamic programming someday so that we
  don't have to do this postprocessing step.
******************************************************************************/
void adjust_starts(struct _gene *genes, int num_genes, struct _node *nodes,
                   int num_nodes, double start_weight)
{
  int i = 0;
  int j = 0;
  /* Tracking variables for indices of top 3 starts */
  int index = 0;
  int tmp_index = 0;
  int max_index[2] = {0};
  /* Score tracking variables */
  double score = 0.0;
  double max_score[2] = {0};
  int valid_distance = 0;
  /* Variables to track modifications to the score for intergenic distance */
  double ig_mod = 0.0;
  double max_ig_mod[2] = {0};

  for (i = 0; i < num_genes; i++)
  {
    /* For each gene, record the index, base score, and modification
       due to intergenic distance with the gene upstream from the start. */
    index = genes[i].start_index;
    score = nodes[index].sscore + nodes[index].cscore;
    ig_mod = get_intergenic_score(genes, num_genes, i, nodes, index,
                                  start_weight, &valid_distance);

    /* Find the best two alternative starts and store their indices, */
    /* scores, and intergenic distance scores in three arrays. */
    get_best_two_alternative_starts(genes, num_genes, i, nodes, num_nodes,
                                    index, start_weight, max_index, max_score,
                                    max_ig_mod);

    /* Dueling starts: See if one start is the clear victor. */

    /* Now we look at the 2nd and 3rd best starts to see if we can find a */
    /* reason to make them our preferred start. */
    /* Change the start if it's a TTG with better coding/RBS/upstream score */
    /* Also change the start if it's <=15bp but has better coding/RBS       */
    for (j = 0; j < 2; j++)
    {
      tmp_index = max_index[j];
      if (tmp_index == -1)
      {
        continue;
      }

      /* Start of less common type but with better coding, rbs, and */
      /* upstream.  Must be 18 or more bases away from original.    */
      if (nodes[tmp_index].subtype != nodes[index].subtype &&
          nodes[tmp_index].tscore < nodes[index].tscore && max_score[j] -
          nodes[tmp_index].tscore >= score-nodes[index].tscore +
          start_weight && nodes[tmp_index].rscore > nodes[index].rscore &&
          nodes[tmp_index].bscore > nodes[index].bscore &&
          nodes[tmp_index].cscore > nodes[index].cscore &&
          abs(nodes[tmp_index].index - nodes[index].index) > 15)
      {
        max_score[j] += nodes[index].tscore-nodes[tmp_index].tscore;
      }

      /* Close starts.  Ignore coding and see if start has better rbs */
      /* and type. */
      else if (abs(nodes[tmp_index].index - nodes[index].index) <= 15 &&
               nodes[tmp_index].rscore + nodes[tmp_index].tscore >
               nodes[index].rscore + nodes[index].tscore &&
               nodes[index].edge == 0 && nodes[tmp_index].edge == 0)
      {
        if (nodes[index].cscore > nodes[tmp_index].cscore)
        {
          max_score[j] += nodes[index].cscore - nodes[tmp_index].cscore;
        }
        if (nodes[index].bscore > nodes[tmp_index].bscore)
        {
          max_score[j] += nodes[index].bscore - nodes[tmp_index].bscore;
        }
      }
      else
      {
        max_score[j] = -1000.0;
      }
    }

    /* Change the gene coordinates to the new maximum. */
    tmp_index = -1;
    for (j = 0; j < 2; j++)
    {
      if (max_index[j] == -1)
      {
        continue;
      }
      if (tmp_index == -1 && max_score[j] + max_ig_mod[j] > score+ig_mod)
      {
        tmp_index = j;
      }
      else if (tmp_index >= 0 && max_score[j] + max_ig_mod[j] >
               max_score[tmp_index] + max_ig_mod[tmp_index])
      {
        tmp_index = j;
      }
    }
    if (tmp_index != -1 && nodes[max_index[tmp_index]].strand == 1)
    {
      nodes[genes[i].start_index].status = 0;
      genes[i].start_index = max_index[tmp_index];
      genes[i].begin = nodes[max_index[tmp_index]].index + 1;
      nodes[max_index[tmp_index]].status = 1;
    }
    else if (tmp_index != -1 && nodes[max_index[tmp_index]].strand == -1)
    {
      nodes[genes[i].start_index].status = 0;
      genes[i].start_index = max_index[tmp_index];
      genes[i].end = nodes[max_index[tmp_index]].index + 1;
      nodes[max_index[tmp_index]].status = 1;
    }
  }
}

/* Given a particular gene and an address of a start node which might */
/* be the start for that gene, return the intergenic distance score */
/* between that start and the gene closest to that start.  Set valid */
/* distance to 0 if this start would cause an unacceptable overlap. */
double get_intergenic_score(struct _gene *genes, int num_genes,
                            int gene_ndx, struct _node *nodes,
                            int node_ndx, double start_weight,
                            int *valid_distance)
{
  *valid_distance = 1;
  if (gene_ndx > 0 && nodes[node_ndx].strand == 1 &&
      nodes[genes[gene_ndx-1].start_index].strand == 1)
  {
    if (nodes[genes[gene_ndx-1].stop_index].index - nodes[node_ndx].index >
        MAX_SAM_OVLP)
    {
      *valid_distance = 0;
      return 0.0;
    }
    else
    {
      return intergenic_mod(&nodes[genes[gene_ndx-1].stop_index],
                            &nodes[node_ndx], start_weight);
    }
  }
  else if (gene_ndx > 0 && nodes[node_ndx].strand == 1 &&
           nodes[genes[gene_ndx-1].start_index].strand == -1)
  {
    if (nodes[genes[gene_ndx-1].start_index].index -
        nodes[node_ndx].index >= 0)
    {
      *valid_distance = 0;
      return 0.0;
    }
    else
    {
      return intergenic_mod(&nodes[genes[gene_ndx-1].start_index],
                            &nodes[node_ndx], start_weight);
    }
  }
  else if (gene_ndx < num_genes-1 && nodes[node_ndx].strand == -1 &&
           nodes[genes[gene_ndx+1].start_index].strand == 1)
  {
    if (nodes[node_ndx].index -
        nodes[genes[gene_ndx+1].start_index].index >= 0)
    {
      *valid_distance = 0;
      return 0.0;
    }
    else
    {
      return intergenic_mod(&nodes[node_ndx],
                            &nodes[genes[gene_ndx+1].start_index],
                            start_weight);
    }
  }
  else if (gene_ndx < num_genes-1 && nodes[node_ndx].strand == -1 &&
           nodes[genes[gene_ndx+1].start_index].strand == -1)
  {
    if (nodes[node_ndx].index -
        nodes[genes[gene_ndx+1].stop_index].index >= 0)
    {
      *valid_distance = 0;
      return 0.0;
    }
    else
    {
      return intergenic_mod(&nodes[node_ndx],
                            &nodes[genes[gene_ndx+1].stop_index],
                            start_weight);
    }
  }
  else
  {
    return 0.0;
  }
}

/* This function finds the best two alternative starts and puts their */
/* indices into the "max_index" array.  If two positive-scoring starts */
/* can't be located, stores a -1 in the index array.  Also stores */
/* scoring information in the max_score and max_ig_mod arrays. */
void get_best_two_alternative_starts(struct _gene *genes, int num_genes,
                                     int gene_ndx, struct _node *nodes,
                                     int num_nodes, int node_ndx,
                                     double start_weight, int *max_index,
                                     double *max_score, double *max_ig_mod)
{
  int valid_distance = 0;
  int i = 0;
  double tmp_ig_mod = 0.0;

  max_index[0] = -1;
  max_index[1] = -1;
  max_score[0] = 0.0;
  max_score[1] = 0.0;
  max_ig_mod[0] = 0.0;
  max_ig_mod[1] = 0.0;

  /* Look up and down 100 nodes to find all the close starts */
  for (i = node_ndx-100; i < node_ndx+100; i++)
  {
    if (i < 0 || i >= num_nodes || i == node_ndx)
    {
      continue;
    }
    if (nodes[i].type == STOP ||
        nodes[i].stop_val != nodes[node_ndx].stop_val)
    {
      continue;
    }
    tmp_ig_mod = get_intergenic_score(genes, num_genes, gene_ndx, nodes, i,
                                      start_weight, &valid_distance);
    if (valid_distance == 0)
    {
      continue;
    }

    /* If it's one of the top two alternatives, add it to our array */
    if (max_index[0] == -1)
    {
      max_index[0] = i;
      max_score[0] = nodes[i].cscore + nodes[i].sscore;
      max_ig_mod[0] = tmp_ig_mod;
    }
    else if (nodes[i].cscore + nodes[i].sscore + tmp_ig_mod > max_score[0])
    {
      max_index[1] = max_index[0];
      max_score[1] = max_score[0];
      max_ig_mod[1] = max_ig_mod[0];
      max_index[0] = i;
      max_score[0] = nodes[i].cscore + nodes[i].sscore;
      max_ig_mod[0] = tmp_ig_mod;
    }
    else if (max_index[1] == -1 || nodes[i].cscore + nodes[i].sscore +
             tmp_ig_mod > max_score[1])
    {
      max_index[1] = i;
      max_score[1] = nodes[i].cscore + nodes[i].sscore;
      max_ig_mod[1] = tmp_ig_mod;
    }
  }
}

/* This function records the metadata associated with a gene into */
/* two text strings: gene_data for start type/RBS motif/etc. and score_data */
/* for the actual numerical scores. */
void record_gene_data(struct _gene *genes, struct _gene_data *gene_data,
                      int num_genes, struct _node *nodes,
                      struct _training *train_data, int seq_counter)
{
  int i = 0;
  int beg_node = 0;
  int end_node = 0;
  int partial_left = 0;
  int partial_right = 0;
  int start_type = 0;
  int stop_type = 0;
  double rbs1 = 0.0;
  double rbs2 = 0.0;
  double confidence = 0.0;
  char sd_string[28][100] = {{0}};
  char sd_spacer[28][20] = {{0}};
  char motif[10] = "";
  char start_string[5][20] = { "ATG", "GTG", "TTG" , "Nonstandard", "Edge" };
  char stop_string[5][20] = { "TAA", "TAG", "TGA" , "Nonstandard", "Edge" };

  /* Initialize RBS string information for default SD */
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

  /* Go through each gene and record its information */
  for (i = 0; i < num_genes; i++)
  {
    beg_node = genes[i].start_index;
    end_node = genes[i].stop_index;

    /* Record basic gene data */
    if ((nodes[beg_node].edge == 1 && nodes[beg_node].strand == 1) ||
        (nodes[end_node].edge == 1 && nodes[beg_node].strand == -1))
    {
      partial_left = 1;
    }
    else
    {
      partial_left = 0;
    }
    if ((nodes[end_node].edge == 1 && nodes[beg_node].strand == 1) ||
        (nodes[beg_node].edge == 1 && nodes[beg_node].strand == -1))
    {
      partial_right = 1;
    }
    else
    {
      partial_right = 0;
    }
    start_type = nodes[beg_node].subtype;
    stop_type = nodes[end_node].subtype;

    sprintf(gene_data[i].gene_data,
            "ID=%d_%d;partial=%d%d;start_type=%s;stop_type=%s;",
            seq_counter, i + 1, partial_left, partial_right,
            start_string[start_type], stop_string[stop_type]);

    /* Record rbs data */
    rbs1 = train_data->rbs_wt[nodes[beg_node].rbs[0]]*train_data->start_weight;
    rbs2 = train_data->rbs_wt[nodes[beg_node].rbs[1]]*train_data->start_weight;
    if (train_data->uses_sd == 1)
    {
      if (rbs1 > rbs2)
      {
        sprintf(gene_data[i].gene_data, "%srbs_motif=%s;rbs_spacer=%s",
                gene_data[i].gene_data, sd_string[nodes[beg_node].rbs[0]],
                sd_spacer[nodes[beg_node].rbs[0]]);
      }
      else
      {
        sprintf(gene_data[i].gene_data, "%srbs_motif=%s;rbs_spacer=%s",
                gene_data[i].gene_data, sd_string[nodes[beg_node].rbs[1]],
                sd_spacer[nodes[beg_node].rbs[1]]);
      }
    }
    else
    {
      mer_text(motif, nodes[beg_node].mot.len, nodes[beg_node].mot.index);
      if (train_data->no_mot > -0.5 && rbs1 > rbs2 &&
          rbs1 > nodes[beg_node].mot.score * train_data->start_weight)
      {
        sprintf(gene_data[i].gene_data, "%srbs_motif=%s;rbs_spacer=%s",
                gene_data[i].gene_data, sd_string[nodes[beg_node].rbs[0]],
                sd_spacer[nodes[beg_node].rbs[0]]);
      }
      else if (train_data->no_mot > -0.5 && rbs2 >= rbs1 &&
               rbs2 > nodes[beg_node].mot.score * train_data->start_weight)
      {
        sprintf(gene_data[i].gene_data, "%srbs_motif=%s;rbs_spacer=%s",
                gene_data[i].gene_data, sd_string[nodes[beg_node].rbs[1]],
                sd_spacer[nodes[beg_node].rbs[1]]);
      }
      else if (nodes[beg_node].mot.len == 0)
      {
        sprintf(gene_data[i].gene_data, "%srbs_motif=None;rbs_spacer=None",
                gene_data[i].gene_data);
      }
      else
      {
        sprintf(gene_data[i].gene_data, "%srbs_motif=%s;rbs_spacer=%dbp",
                gene_data[i].gene_data, motif, nodes[beg_node].mot.spacer);
      }
    }
    sprintf(gene_data[i].gene_data, "%s;gc_cont=%.3f", gene_data[i].gene_data,
            nodes[beg_node].gc_cont);

    /* Record score data */
    confidence = calculate_confidence(nodes[beg_node].cscore +
                                      nodes[beg_node].sscore,
                                      train_data->start_weight);
    sprintf(gene_data[i].score_data,
      "conf=%.2f;score=%.2f;cscore=%.2f;sscore=%.2f;rscore=%.2f;bscore=%.2f;",
      confidence, nodes[beg_node].cscore+nodes[beg_node].sscore,
      nodes[beg_node].cscore, nodes[beg_node].sscore, nodes[beg_node].rscore,
      nodes[beg_node].bscore);
    sprintf(gene_data[i].score_data, "%stscore=%.2f",
            gene_data[i].score_data, nodes[beg_node].tscore);
  }
}

/* Print the genes.  The 'format' variable indicates which format to use. */
void print_genes(FILE *fp, struct _gene *genes, struct _gene_data *gene_data,
                 int num_genes, struct _node *nodes, int seq_length,
                 int format, int seq_counter, int mode, char *preset_data,
                 struct _training *train_data, char *header, char *short_hdr,
                 char *version)
{
  int i = 0;
  int beg_node = 0;
  int end_node = 0;
  char left[50] = "";
  char right[50] = "";
  char seq_data[MAX_LINE*2] = "";
  char run_data[MAX_LINE] = "";

  /* Initialize sequence data */
  sprintf(seq_data, "seqnum=%d;seqlen=%d;seqhdr=\"%s\"", seq_counter,
          seq_length, header);

  /* Initialize run data string */
  if (mode == MODE_ANON)
  {
    sprintf(run_data, "version=Prodigal.v%s;run_type=Anonymous;", version);
    sprintf(run_data, "%smodel=\"%s\";", run_data, preset_data);
  }
  else
  {
    sprintf(run_data, "version=Prodigal.v%s;run_type=Normal;", version);
    sprintf(run_data, "%smodel=\"Ab initio\";", run_data);
  }
  sprintf(run_data, "%sgc_cont=%.2f;transl_table=%d;uses_sd=%d", run_data,
          train_data->gc*100.0, train_data->trans_table, train_data->uses_sd);

  strcpy(left, "");
  strcpy(right, "");

  /* Print the gff header once */
  if (format == 3 && seq_counter == 1)
  {
    fprintf(fp, "##gff-version  3\n");
  }

  /* Print sequence/model information */
  if (format == 1)
  {
    fprintf(fp, "DEFINITION  %s;%s\n", seq_data, run_data);
    fprintf(fp, "FEATURES             Location/Qualifiers\n");
  }
  else if (format > 1 && format < 4)
  {
    fprintf(fp, "# Sequence Data: %s\n", seq_data);
    fprintf(fp, "# Model Data: %s\n", run_data);
  }
  else if (format == 4)
  {
    fprintf(fp, ">Feature %s\n", short_hdr);
    fprintf(fp, "1\t%d\tREFERENCE\n", seq_length);
  }

  /* Print the genes */
  for (i = 0; i < num_genes; i++)
  {
    beg_node = genes[i].start_index;
    end_node = genes[i].stop_index;

    /* Print the coordinates and data */
    if (nodes[beg_node].strand == 1)
    {
      if (nodes[beg_node].edge == 1)
      {
        sprintf(left, "<%d", genes[i].begin);
      }
      else
      {
        sprintf(left, "%d", genes[i].begin);
      }
      if (nodes[end_node].edge == 1)
      {
        sprintf(right, ">%d", genes[i].end);
      }
      else
      {
        sprintf(right, "%d", genes[i].end);
      }

      if (format == 1)
      {
        fprintf(fp, "     CDS             %s..%s\n", left, right);
        fprintf(fp, "                     ");
        fprintf(fp, "/note=\"%s;%s\"\n", gene_data[i].gene_data,
                gene_data[i].score_data);
      }
      if (format == 2) fprintf(fp, ">%d_%d_%d_+\n", i+1, genes[i].begin,
          genes[i].end);
      {
      }
      if (format == 3)
      {
        fprintf(fp, "%s\tProdigal_v%s\tCDS\t%d\t%d\t%.1f\t+\t0\t%s;%s",
                short_hdr, version, genes[i].begin, genes[i].end,
                nodes[beg_node].cscore + nodes[beg_node].sscore,
                gene_data[i].gene_data, gene_data[i].score_data);
        fprintf(fp, "\n");
      }
      if (format == 4)
      {
        fprintf(fp, "%s\t%s\tCDS\n", left, right);
        if (nodes[beg_node].edge == 1)
        {
          fprintf(fp, "\t\t\tcodon_start\t%d\n", genes[i].begin);
        }
        fprintf(fp, "\t\t\tnote\t%s;%s\n", gene_data[i].gene_data,
                gene_data[i].score_data);
      }
    }
    else
    {
      if (nodes[end_node].edge == 1)
      {
        sprintf(left, "<%d", genes[i].begin);
      }
      else
      {
        sprintf(left, "%d", genes[i].begin);
      }
      if (nodes[beg_node].edge == 1)
      {
        sprintf(right, ">%d", genes[i].end);
      }
      else
      {
        sprintf(right, "%d", genes[i].end);
      }

      if (format == 1)
      {
        fprintf(fp, "     CDS             complement(%s..%s)\n", left, right);
        fprintf(fp, "                     ");
        fprintf(fp, "/note=\"%s;%s\"\n", gene_data[i].gene_data,
                gene_data[i].score_data);
      }
      if (format == 2)
      {
        fprintf(fp, ">%d_%d_%d_-\n", i+1, genes[i].begin, genes[i].end);
      }
      if (format == 3)
      {
        fprintf(fp, "%s\tProdigal_v%s\tCDS\t%d\t%d\t%.1f\t-\t0\t%s;%s",
                short_hdr, version, genes[i].begin, genes[i].end,
                nodes[beg_node].cscore + nodes[beg_node].sscore,
                gene_data[i].gene_data, gene_data[i].score_data);
        fprintf(fp, "\n");
      }
      if (format == 4)
      {
        fprintf(fp, "%s\t%s\tCDS\n", right, left);
        if (nodes[beg_node].edge == 1)
        {
          fprintf(fp, "\t\t\tcodon_start\t%d\n", genes[i].end);
        }
        fprintf(fp, "\t\t\tnote\t%s;%s\n", gene_data[i].gene_data,
                gene_data[i].score_data);
      }
    }
  }

  /* Genbank footer */
  if (format == 1)
  {
    fprintf(fp, "//\n");
  }
}

/* Print the gene translations */
void write_translations(FILE *fh, struct _gene *genes,
                        struct _gene_data *gene_data, int num_genes,
                        struct _node *nodes, unsigned char *seq,
                        unsigned char *rseq, unsigned char *useq,
                        int seq_length, int trans_table, char *short_hdr)
{
  int i = 0;
  int j = 0;

  if (fh == NULL)
  {
    return;
  }

  for (i = 0; i < num_genes; i++)
  {
    if (nodes[genes[i].start_index].strand == 1)
    {
      fprintf(fh, ">%s_%d # %d # %d # 1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, gene_data[i].gene_data);
      for (j = genes[i].begin; j < genes[i].end; j+=3)
      {
        if (is_n(useq, j-1) == 1 || is_n(useq, j) == 1 || is_n(useq, j+1) == 1)
        {
          fprintf(fh, "X");
        }
        else
        {
          fprintf(fh, "%c", amino(seq, j-1, trans_table,
                  (j==genes[i].begin?1:0) &&
                  (1 - nodes[genes[i].start_index].edge)));
        }
        if ((j - genes[i].begin) % 180 == 177)
        {
          fprintf(fh, "\n");
        }
      }
      if ((j - genes[i].begin) % 180 != 0)
      {
        fprintf(fh, "\n");
      }
    }
    else
    {
      fprintf(fh, ">%s_%d # %d # %d # -1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, gene_data[i].gene_data);
      for (j = seq_length + 1 - genes[i].end; j < seq_length + 1 -
           genes[i].begin; j+=3)
      {
        if (is_n(useq, seq_length - j) == 1 || is_n(useq, seq_length - 1 - j)
            == 1 || is_n(useq, seq_length - 2 - j) == 1)
        {
          fprintf(fh, "X");
        }
        else
        {
          fprintf(fh, "%c", amino(rseq, j-1, trans_table,
                  (j==seq_length + 1 - genes[i].end?1:0) &&
                  (1 - nodes[genes[i].start_index].edge)));
        }
        if ((j - seq_length - 1 + genes[i].end) % 180 == 177)
        {
          fprintf(fh, "\n");
        }
      }
      if ((j - seq_length - 1 + genes[i].end) % 180 != 0)
      {
        fprintf(fh, "\n");
      }
    }
  }
}

/* Print the gene nucleotide sequences */
void write_nucleotide_seqs(FILE *fh, struct _gene *genes,
                           struct _gene_data *gene_data, int num_genes,
                           struct _node *nodes, unsigned char *seq,
                           unsigned char *rseq, unsigned char *useq,
                           int seq_length, char *short_hdr)
{
  int i = 0;
  int j = 0;

  if (fh == NULL)
  {
    return;
  }

  for (i = 0; i < num_genes; i++)
  {
    if (nodes[genes[i].start_index].strand == 1)
    {
      fprintf(fh, ">%s_%d # %d # %d # 1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, gene_data[i].gene_data);
      for (j = genes[i].begin-1; j < genes[i].end; j++)
      {
        if (is_n(useq, j) == 1)
        {
          fprintf(fh, "N");
        }
        else if (is_a(seq, j) == 1)
        {
          fprintf(fh, "A");
        }
        else if (is_t(seq, j) == 1)
        {
          fprintf(fh, "T");
        }
        else if (is_g(seq, j) == 1)
        {
          fprintf(fh, "G");
        }
        else if (is_c(seq, j) == 1)
        {
          fprintf(fh, "C");
        }
        if ((j - genes[i].begin + 1) % 70 == 69)
        {
          fprintf(fh, "\n");
        }
      }
      if ((j - genes[i].begin + 1) % 70 != 0)
      {
        fprintf(fh, "\n");
      }
    }
    else
    {
      fprintf(fh, ">%s_%d # %d # %d # -1 # %s\n", short_hdr, i+1,
              genes[i].begin, genes[i].end, gene_data[i].gene_data);
      for (j = seq_length-genes[i].end; j < seq_length+1-genes[i].begin; j++)
      {
        if (is_n(useq, seq_length - 1 - j) == 1)
        {
          fprintf(fh, "N");
        }
        else if (is_a(rseq, j) == 1)
        {
          fprintf(fh, "A");
        }
        else if (is_t(rseq, j) == 1)
        {
          fprintf(fh, "T");
        }
        else if (is_g(rseq, j) == 1)
        {
          fprintf(fh, "G");
        }
        else if (is_c(rseq, j) == 1)
        {
          fprintf(fh, "C");
        }
        if ((j - seq_length + genes[i].end) % 70 == 69)
        {
          fprintf(fh, "\n");
        }
      }
      if ((j - seq_length + genes[i].end) % 70 != 0)
      {
        fprintf(fh, "\n");
      }
    }
  }
}

/******************************************************************************
  Write detailed scoring information about every single possible gene.  Only
  done at the user's request.
******************************************************************************/
void write_start_file(FILE *fh, struct _node *nodes, int num_nodes,
                      struct _training *train_data, int seq_counter,
                      int seq_length, int mode, char *preset_data,
                      char *version, char *header)
{
  int i = 0;
  int prev_stop = -1;
  int prev_strand = 0;
  int start_type = 0;
  double rbs1 = 0.0;
  double rbs2 = 0.0;
  char sd_string[28][100] = {{0}};
  char sd_spacer[28][20] = {{0}};
  char motif[10] = "";
  char start_string[5][20] = { "ATG", "GTG", "TTG" , "Nonstandard", "Edge" };
  char seq_data[MAX_LINE*2] = "";
  char run_data[MAX_LINE] = "";

  if (fh == NULL)
  {
    return;
  }

  /* Initialize sequence data */
  sprintf(seq_data, "seqnum=%d;seqlen=%d;seqhdr=\"%s\"", seq_counter,
          seq_length, header);

  /* Initialize run data string */
  if (mode == MODE_ANON)
  {
    sprintf(run_data, "version=Prodigal.v%s;run_type=Anonymous;", version);
    sprintf(run_data, "%smodel=\"%s\";", run_data, preset_data);
  }
  else
  {
    sprintf(run_data, "version=Prodigal.v%s;run_type=Normal;", version);
    sprintf(run_data, "%smodel=\"Ab initio\";", run_data);
  }
  sprintf(run_data, "%sgc_cont=%.2f;transl_table=%d;uses_sd=%d", run_data,
          train_data->gc*100.0, train_data->trans_table, train_data->uses_sd);

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

  qsort(nodes, num_nodes, sizeof(struct _node), &stopcmp_nodes);

  fprintf(fh, "# Sequence Data: %s\n", seq_data);
  fprintf(fh, "# Run Data: %s\n\n", run_data);

  fprintf(fh, "Beg\tEnd\tStd\tTotal\tCodPot\tLenSc\tStrtSc\tCodon\tRBSMot\t");
  fprintf(fh, "Spacer\tRBSScr\tUpsScr\tTypeScr\tStpScr\tGCCont\n");
  for (i = 0; i < num_nodes; i++)
  {
    if (nodes[i].type == STOP)
    {
      continue;
    }
    start_type = nodes[i].subtype;
    if (nodes[i].stop_val != prev_stop || nodes[i].strand != prev_strand)
    {
      prev_stop = nodes[i].stop_val;
      prev_strand = nodes[i].strand;
      fprintf(fh, "\n");
    }
    if (nodes[i].strand == 1)
    {
      fprintf(fh, "%d\t%d\t+\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t", nodes[i].index+1,
              nodes[i].stop_val+3, nodes[i].cscore+nodes[i].sscore,
              nodes[i].pscore, nodes[i].lscore, nodes[i].sscore,
              start_string[start_type]);
    }
    if (nodes[i].strand == -1)
    {
      fprintf(fh, "%d\t%d\t-\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t",
              nodes[i].stop_val-1, nodes[i].index+1,
              nodes[i].cscore+nodes[i].sscore, nodes[i].pscore,
              nodes[i].lscore, nodes[i].sscore, start_string[start_type]);
    }
    rbs1 = train_data->rbs_wt[nodes[i].rbs[0]]*train_data->start_weight;
    rbs2 = train_data->rbs_wt[nodes[i].rbs[1]]*train_data->start_weight;
    if (train_data->uses_sd == 1)
    {
      if (rbs1 > rbs2)
      {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nodes[i].rbs[0]],
                sd_spacer[nodes[i].rbs[0]], nodes[i].rscore);
      }
      else
      {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nodes[i].rbs[1]],
                sd_spacer[nodes[i].rbs[1]], nodes[i].rscore);
      }
    }
    else
    {
      mer_text(motif, nodes[i].mot.len, nodes[i].mot.index);
      if (train_data->no_mot > -0.5 && rbs1 > rbs2 &&
          rbs1 > nodes[i].mot.score * train_data->start_weight)
      {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nodes[i].rbs[0]],
                sd_spacer[nodes[i].rbs[0]], nodes[i].rscore);
      }
      else if (train_data->no_mot > -0.5 && rbs2 >= rbs1 &&
               rbs2 > nodes[i].mot.score * train_data->start_weight)
      {
        fprintf(fh, "%s\t%s\t%.2f\t", sd_string[nodes[i].rbs[1]],
                sd_spacer[nodes[i].rbs[1]], nodes[i].rscore);
      }
      else
      {
        if (nodes[i].mot.len == 0)
        {
          fprintf(fh, "None\tNone\t%.2f\t", nodes[i].rscore);
        }
        else
        {
          fprintf(fh, "%s\t%dbp\t%.2f\t", motif, nodes[i].mot.spacer,
                  nodes[i].rscore);
        }
      }
    }
    fprintf(fh, "%.2f\t%.2f\t%.3f\n", nodes[i].bscore, nodes[i].tscore,
            nodes[i].gc_cont);
  }
  fprintf(fh, "\n");
  qsort(nodes, num_nodes, sizeof(struct _node), &compare_nodes);
}

/* Convert score to a percent confidence */
double calculate_confidence(double score, double start_weight)
{
  double conf = 0.0;

  if (score/start_weight < 41)
  {
    conf = exp(score/start_weight);
    conf = (conf/(conf+1)) * 100.0;
  }
  else
  {
    conf = 100.00;
  }
  if (conf <= 50.00)
  {
    conf = 50.00;
  }
  return conf;
}

/******************************************************************************
  Make node status flag match what we have in the gene list.  Necessary
  mainly for anonymous runs where we don't have access to the dynamic
  programming pointers.
******************************************************************************/
void match_nodes_to_genes(struct _gene *genes, int num_genes,
                          struct _node *nodes, int num_nodes)
{
  int i = 0;

  /* First eliminate all nodes */
  for (i = 0; i < num_nodes; i++)
  {
    nodes[i].status = 0;
  }

  /* Second, mark all the nodes in the gene list as status 1 */
  for (i = 0; i < num_genes; i++)
  {
    nodes[genes[i].start_index].status = 1;
    nodes[genes[i].stop_index].status = 1;
  }
}
