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

#include "setup.h"
#include "summary.h"

int main(int argc, char *argv[])
{

  int slen, nn, ng, i, ipath, outfmt, max_phase;
  int closed, num_seq, cross_gaps, quiet;
  int piped, max_slen, genetic_code, gene_len_flag;
  int mode; /* 0 = normal, 1 = training, 2 = anonymous/metagenomic */
  int force_nonsd = 0; /* deprecated and slated for removal PDH */
  double max_score, gc, low, high;
  unsigned char *seq, *rseq, *useq;
  char train_file[MAX_LINE], start_file[MAX_LINE], amino_file[MAX_LINE];
  char nuc_file[MAX_LINE], summ_file[MAX_LINE], input_file[MAX_LINE];
  char output_file[MAX_LINE];
  char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
  FILE *input_ptr, *output_ptr, *start_ptr, *amino_ptr, *nuc_ptr, *summ_ptr;
  struct _node *nodes;
  struct _gene *genes, **anon_genes;
  struct _gene_data *gene_data;
  struct _training tinf = {0};
  struct _preset_genome_bin presets[NUM_PRESET_GENOME] = {{0}};
  struct _summary statistics = {0};

  /* Allocate memory for data structures */
  if (initialize_data_structures(&seq, &rseq, &useq, &nodes, &genes,
                                 &gene_data, presets, &anon_genes) == -1)
  {
    perror("\nError: Dynamic memory allocation failed.");
    exit(1);
  }

  /* Initialize all variables to 0 and pointers to NULL */
  nn = 0;
  slen = 0;
  ipath = 0;
  ng = 0;
  num_seq = 0;
  max_phase = 0;
  max_score = -100.0;
  piped = 0;
  max_slen = 0;
  gene_len_flag = 0;
  input_ptr = stdin;
  output_ptr = stdout;
  start_ptr = NULL;
  amino_ptr = NULL;
  nuc_ptr = NULL;
  summ_ptr = NULL;

  /* Parse and validate the command line arguments */
  parse_arguments(argc, argv, input_file, output_file, train_file,
                  amino_file, nuc_file, start_file, summ_file, &mode,
                  &outfmt, &genetic_code, &closed, &cross_gaps, &quiet);

  /* Defaults for genetic code, output format, start weight */
  tinf.start_weight = 4.35;
  if (genetic_code == -1)
  {
    tinf.trans_table = 11; /* 11 default genetic code */
  }
  else
  {
    tinf.trans_table = genetic_code;
  }
  if (outfmt == -1)
  {
    outfmt = 3; /* GFF default output format */
  }

  /* Print header */
  if (quiet == 0)
  {
    header(mode);
  }

  /* Look for input on stdin and handle Windows' inability to rewind stdin */
  if (mode == 0 && train_file[0] == '\0' && input_file[0] == '\0')
  {
    piped = detect_input_and_handle_windows_stdin(argc, quiet, input_file);
  }

  /* Read in the training file (if specified) */
  if (train_file[0] != '\0')
  {
    if (quiet == 0)
    {
      fprintf(stderr, "Reading in training data from file %s...", train_file);
    }
    if (read_training_file(train_file, &tinf) == -1)
    {
      perror("\n\nError: training file did not read correctly!");
      exit(6);
    }
    if (quiet == 0)
    {
      fprintf(stderr, "done.\n");
      fprintf(stderr, "-------------------------------------\n");
    }
  }

  /* Check i/o files (if specified) and prepare them for reading/writing */
  open_files(input_file, output_file, start_file, amino_file, nuc_file,
             summ_file, &input_ptr, &output_ptr, &start_ptr, &amino_ptr,
             &nuc_ptr, &summ_ptr);

  /***************************************************************************
    Single Genome Training:  Read in the sequence(s) and perform the
    training on them.
  ***************************************************************************/
  if (mode == 1 || (mode == 0 && train_file[0] == '\0'))
  {
    if (quiet == 0)
    {
      fprintf(stderr, "Reading in the sequence(s) to train...");
    }
    slen = read_seq_training(input_ptr, seq, useq, &(tinf.gc), closed,
                             &statistics.num_contig);
    rcom_seq(seq, rseq, useq, slen);
    calc_avg_train_contig_len(&statistics, slen);
    if (quiet == 0)
    {
      fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen,
              tinf.gc * 100.0);
    }

    /* Grab more memory if sequence is larger than our default allocation */
    if (slen > STT_NOD*8)
    {
      nodes = (struct _node *)realloc(nodes,
                                      (int)(slen / 8) * sizeof(struct _node));
      if (nodes == NULL)
      {
        perror("Realloc failed on nodes\n\n");
        exit(11);
      }
    }

    /* Build the training set and score the coding of every start-stop pair */
    if (quiet == 0)
    {
      fprintf(stderr, "Building training set using genetic code %d...",
              tinf.trans_table);
    }
    build_training_set(nodes, &tinf, &statistics, seq, rseq, useq, slen, &nn,
                       closed, cross_gaps);
    if (quiet == 0)
    {
      fprintf(stderr, "done!\n");
    }

    /***********************************************************************
      Check average gene length to see if translation table looks good or
      if there is substantial gene decay.  If no genetic code was specified,
      try genetic code 4 to see if it solves the problem.
    ***********************************************************************/
    if (quiet == 0)
    {
      fprintf(stderr, "Checking average training gene length...");
    }
    gene_len_flag = bad_train_gene_length(statistics);
    if (gene_len_flag == 0)
    { /* ok */
      if (quiet == 0)
      {
        fprintf(stderr, "%.1f, looks ok.\n", statistics.avg_comp_gene_len);
      }
    }
    else if (gene_len_flag == 1)
    { /* too many partial genes */
      if (quiet == 0)
      {
        fprintf(stderr, "low but sequence is really drafty.\n");
      }
      low_gene_len_warning(gene_len_flag, statistics);
    }
    else if (gene_len_flag == 2)
    { /* poor avg gene length */
      if (quiet == 0)
      {
        fprintf(stderr, "%.1f, too low.\n", statistics.avg_comp_gene_len);
      }
      if (genetic_code == -1)
      {
        if (quiet == 0)
        {
          fprintf(stderr, "Trying genetic code 4...");
        }
        tinf.trans_table = 4;
        build_training_set(nodes, &tinf, &statistics, seq, rseq, useq, slen,
                           &nn, closed, cross_gaps);
        gene_len_flag = bad_train_gene_length(statistics);
        if (gene_len_flag < 2)
        {
          if (quiet == 0)
          {
            fprintf(stderr, "looks good, going with genetic code 4.\n");
          }
        }
        else
        {
          if (quiet == 0)
          {
            fprintf(stderr, "still bad, switching back to genetic code 11.\n");
            fprintf(stderr, "Redoing genome with genetic code 11...");
          }
          tinf.trans_table = 11;
          build_training_set(nodes, &tinf, &statistics, seq, rseq, useq, slen,
                             &nn, closed, cross_gaps);
          if (quiet == 0)
          {
            fprintf(stderr, "done.\n");
          }
          low_gene_len_warning(gene_len_flag, statistics);
        }
      }
      else
      {
        low_gene_len_warning(gene_len_flag, statistics);
      }
    }

    /***********************************************************************
      Determine if this organism uses Shine-Dalgarno or not and score the
      nodes appropriately.
    ***********************************************************************/
    if (quiet == 0)
    {
      fprintf(stderr, "Examining upstream regions and training starts...");
    }
    rbs_score(seq, rseq, slen, nodes, nn, tinf.rbs_wt);
    train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
    determine_sd_usage(&tinf);
    if (force_nonsd == 1)
    {
      tinf.uses_sd = 0;
    }
    if (tinf.uses_sd == 0)
    {
      train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf);
    }
    if (quiet == 0)
    {
      fprintf(stderr, "done.\n");
    }

    /* If training specified, write the training file and exit. */
    if (mode == 1)
    {
      if (quiet == 0)
      {
        fprintf(stderr, "Writing data to training file %s...", output_file);
      }
      if (write_training_file(output_file, &tinf) != 0)
      {
        perror("\nError: could not write training file!");
        exit(12);
      }
      else
      {
        if (quiet == 0)
        {
          fprintf(stderr, "done.\n");
        }
        exit(0);
      }
    }

    /* Rewind input file */
    if (quiet == 0)
    {
      fprintf(stderr, "-------------------------------------\n");
    }
    if (fseek(input_ptr, 0, SEEK_SET) == -1)
    {
      perror("\nError: could not rewind input file.");
      exit(13);
    }
    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, (slen /4 + 1) * sizeof(unsigned char));
    memset(rseq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(useq, 0, (slen / 8 + 1) * sizeof(unsigned char));
    memset(nodes, 0, nn * sizeof(struct _node));
    memset(&statistics, 0, sizeof(struct _summary));
    nn = 0;
    slen = 0;
    ipath = 0;
    num_seq = 0;
  }

  /* Initialize the training files for an anonymous request */
  else if (mode == 2)
  {
    if (quiet == 0)
    {
      fprintf(stderr, "Initializing preset training files...");
    }
    initialize_preset_genome_bins(presets);
    if (quiet == 0)
    {
      fprintf(stderr, "done.\n");
      fprintf(stderr, "-------------------------------------\n");
    }
  }

  /* Print out header for gene finding phase */
  if (quiet == 0)
  {
    if (mode == 2)
    {
      fprintf(stderr, "Mode: Anonymous, Phase: Gene Finding\n");
    }
    else
    {
      fprintf(stderr, "Mode: Normal, Phase: Gene Finding\n");
    }
  }

  /* Read and process each sequence in the file in succession */
  sprintf(cur_header, "Prodigal_Seq_1");
  sprintf(new_header, "Prodigal_Seq_2");
  while ((slen = next_seq_multi(input_ptr, seq, useq, &num_seq, &gc,
         cur_header, new_header)) != -1)
  {
    rcom_seq(seq, rseq, useq, slen);

    if (quiet == 0)
    {
      fprintf(stderr, "Finding genes in sequence #%d (%d bp)...", num_seq,
              slen);
    }

    /* Reallocate memory if this is the biggest sequence we've seen */
    if (slen > max_slen && slen > STT_NOD*8)
    {
      nodes = (struct _node *)realloc(nodes,
                                      (int)(slen / 8) * sizeof(struct _node));
      if (nodes == NULL)
      {
        fprintf(stderr, "Realloc failed on nodes\n\n");
        exit(11);
      }
      max_slen = slen;
    }

    /* Calculate short header for this sequence */
    calc_short_header(cur_header, short_header, num_seq);

    if (mode != 2) /* Single Genome Version */
    {

      /***********************************************************************
        Find all the potential starts and stops, sort them, and create a
        comprehensive list of nodes for dynamic programming.
      ***********************************************************************/
      nn = add_nodes(seq, rseq, useq, slen, nodes, closed, cross_gaps,
                     tinf.trans_table);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

      /***********************************************************************
        Second dynamic programming, using the dicodon statistics as the
        scoring function.
      ***********************************************************************/
      score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, mode);
      record_overlapping_starts(nodes, nn, tinf.start_weight, 1);
      ipath = dynamic_programming(nodes, nn, tinf.bias, tinf.start_weight, 1);
      if (start_ptr != NULL)
      {
        write_start_file(start_ptr, nodes, nn, &tinf, num_seq, slen, mode,
                         NULL, VERSION, cur_header);
      }
      eliminate_bad_genes(nodes, ipath, tinf.start_weight);
      ng = add_genes(genes, nodes, ipath);
      adjust_close_starts(genes, ng, nodes, nn, tinf.start_weight);
      record_gene_data(genes, gene_data, ng, nodes, &tinf, num_seq);
      if (quiet == 0)
      {
        fprintf(stderr, "done.\n");
      }

      /* Output the genes */
      print_genes(output_ptr, genes, gene_data, ng, nodes, slen, outfmt,
                  num_seq, mode, NULL, &tinf, cur_header, short_header,
                  VERSION);
      fflush(output_ptr);
      if (amino_ptr != NULL)
      {
        write_translations(amino_ptr, genes, gene_data, ng, nodes, seq, rseq,
                           useq, slen, tinf.trans_table, num_seq,
                           short_header);
      }
      if (nuc_ptr != NULL)
      {
        write_nucleotide_seqs(nuc_ptr, genes, gene_data, ng, nodes, seq, rseq,
                              useq, slen, num_seq, short_header);
      }
    }

    else /* Anonymous (Metagenomic) Version */
    {
      low = 0.88495 * gc - 0.0102337;
      if (low > 0.65)
      {
        low = 0.65;
      }
      high = 0.86596 * gc + .1131991;
      if (high < 0.35)
      {
        high = 0.35;
      }

      max_score = -100.0;
      for (i = 0; i < NUM_PRESET_GENOME; i++)
      {
        if (i == 0 || presets[i].tinf->trans_table !=
            presets[i-1].tinf->trans_table)
        {
          memset(nodes, 0, nn * sizeof(struct _node));
          nn = add_nodes(seq, rseq, useq, slen, nodes, closed, cross_gaps,
                         presets[i].tinf->trans_table);
          qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        }
        if (presets[i].tinf->gc < low || presets[i].tinf->gc > high)
        {
          continue;
        }
        reset_node_scores(nodes, nn);
        score_nodes(seq, rseq, slen, nodes, nn, presets[i].tinf, closed, mode);
        record_overlapping_starts(nodes, nn, presets[i].tinf->start_weight, 1);
        ipath = dynamic_programming(nodes, nn, presets[i].tinf->bias,
                                    presets[i].tinf->start_weight, 1);
        if (nodes[ipath].score > max_score)
        {
          max_phase = i;
          max_score = nodes[ipath].score;
          eliminate_bad_genes(nodes, ipath, presets[i].tinf->start_weight);
          ng = add_genes(anon_genes[i], nodes, ipath);
          adjust_close_starts(anon_genes[i], ng, nodes, nn,
                              presets[i].tinf->start_weight);
          record_gene_data(anon_genes[i], gene_data, ng, nodes,
                           presets[i].tinf, num_seq);
        }
      }

      /* Recover the nodes for the best of the runs */
      memset(nodes, 0, nn * sizeof(struct _node));
      nn = add_nodes(seq, rseq, useq, slen, nodes, closed, cross_gaps,
                     presets[max_phase].tinf->trans_table);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
      score_nodes(seq, rseq, slen, nodes, nn, presets[max_phase].tinf, closed,
                  mode);
      if (start_ptr != NULL)
      {
        write_start_file(start_ptr, nodes, nn, presets[max_phase].tinf,
                         num_seq, slen, mode, presets[max_phase].desc, VERSION,
                         cur_header);
      }

      if (quiet == 0)
      {
        fprintf(stderr, "done.\n");
      }

      /* Output the genes */
      print_genes(output_ptr, anon_genes[max_phase], gene_data, ng, nodes,
                  slen, outfmt, num_seq, mode, presets[max_phase].desc,
                  presets[max_phase].tinf, cur_header, short_header, VERSION);
      fflush(output_ptr);
      if (amino_ptr != NULL)
      {
        write_translations(amino_ptr, anon_genes[max_phase], gene_data, ng,
                           nodes, seq, rseq, useq, slen,
                           presets[max_phase].tinf->trans_table, num_seq,
                           short_header);
      }
      if (nuc_ptr != NULL)
      {
        write_nucleotide_seqs(nuc_ptr, anon_genes[max_phase], gene_data, ng,
                              nodes, seq, rseq, useq, slen, num_seq,
                              short_header);
      }
    }

    /* Reset all the sequence/dynamic programming variables */
    zero_sequence(seq, rseq, useq, slen);
    zero_nodes(nodes, nn);
    nn = 0;
    slen = 0;
    ipath = 0;
    strcpy(cur_header, new_header);
    sprintf(new_header, "Prodigal_Seq_%d\n", num_seq + 1);
  }

  /* Flag an error if we saw no sequences */
  if (num_seq == 0)
  {
    fprintf(stderr, "\nError:  no input sequences to analyze.\n\n");
    exit(17);
  }

  /* Cleanup: free variables, close filehandles, remove tmp file */
  free_variables(seq, rseq, useq, nodes, genes, presets, anon_genes);
  close_filehandles(input_ptr, output_ptr, start_ptr, amino_ptr, nuc_ptr,
                    summ_ptr);
  if (piped == 1 && remove(input_file) != 0)
  {
    fprintf(stderr, "Could not delete tmp file %s.\n", input_file);
    exit(18);
  }

  /* Exit successfully */
  exit(0);
}
