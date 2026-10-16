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

#include <sys/stat.h>
#include <unistd.h>
#include "prodigal_internal.h"
#include "fptr.h"


#define VERSION PRODIGAL_VERSION_STRING
#define DATE "April, 2026"

#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000


#ifndef PRODIGAL_NO_MAIN

void version();
void usage(char *);
void help();
int copy_standard_input_to_file(char *, int);

int main(int argc, char *argv[]) {

  int rv, nn, ng, i, ipath, do_training, output, max_phase;
  int user_tt, num_seq, quiet;
  int piped, fnum;
  double max_score, path_score, gc, low, high;
  char *train_file, *start_file, *trans_file, *nuc_file;
  char *input_file, *output_file, input_copy[MAX_LINE];
  char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
  FILE *output_ptr, *start_ptr, *trans_ptr, *nuc_ptr;
  fptr input_ptr = NULL;
  struct stat fbuf;
  pid_t pid;

  /* Library context and config */
  prodigal_config_t config;
  prodigal_ctx_t *ctx;

  /* Initialize config with defaults */
  prodigal_config_init(&config);

  nn = 0; ipath = 0; ng = 0;
  user_tt = 0; num_seq = 0; quiet = 0;
  max_phase = 0; max_score = -100.0;
  train_file = NULL; do_training = 0;
  start_file = NULL; trans_file = NULL; nuc_file = NULL;
  start_ptr = stdout; trans_ptr = stdout; nuc_ptr = stdout;
  input_file = NULL; output_file = NULL; piped = 0;
  output_ptr = stdout;
  output = 0;

  /* Filename for input copy if needed */
  pid = getpid();
  sprintf(input_copy, "tmp.prodigal.stdin.%d", pid);

  /* Parse the command line arguments */
  for(i = 1; i < argc; i++) {
    if(i == argc-1 && (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "-T") == 0
       || strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0 ||
       strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "-g") == 0 ||
       strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-F") == 0 ||
       strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0 ||
       strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0 ||
       strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-O") == 0 ||
       strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-P") == 0))
      usage("-a/-f/-g/-i/-o/-p/-s options require parameters.");
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "-C") == 0)
      config.closed_ends = 1;
    else if(strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "-Q") == 0)
      quiet = 1;
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-M") == 0)
      config.mask_regions = 1;
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "-N") == 0)
      config.force_nonsd = 1;
    else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-H") == 0) help();
    else if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "-V") == 0) version();
    else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0) {
      trans_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "-d") == 0) {
      nuc_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0) {
      input_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-O") == 0) {
      output_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0) {
      start_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "-T") == 0) {
      train_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "-G") == 0) {
      config.trans_table = atoi(argv[i+1]);
      if(config.trans_table < 1 || config.trans_table > 25 ||
         config.trans_table == 7 || config.trans_table == 8 ||
         (config.trans_table >= 17 && config.trans_table <= 20))
        usage("Invalid translation table specified.");
      user_tt = config.trans_table;
      i++;
    }
    else if(strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-P") == 0) {
      if(argv[i+1][0] == '0' || argv[i+1][0] == 's' || argv[i+1][0] ==
              'S') config.meta_mode = 0;
      else if(argv[i+1][0] == '1' || argv[i+1][0] == 'm' || argv[i+1][0] ==
              'M') config.meta_mode = 1;
      else usage("Invalid meta/single genome type specified.");
      i++;
    }
    else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-F") == 0) {
      if(strncmp(argv[i+1], "0", 1) == 0 || strcmp(argv[i+1], "gbk") == 0 ||
         strcmp(argv[i+1], "GBK") == 0)
        output = 0;
      else if(strncmp(argv[i+1], "1", 1) == 0 || strcmp(argv[i+1], "gca") == 0
              || strcmp(argv[i+1], "GCA") == 0)
        output = 1;
      else if(strncmp(argv[i+1], "2", 1) == 0 || strcmp(argv[i+1], "sco") == 0
              || strcmp(argv[i+1], "SCO") == 0)
        output = 2;
      else if(strncmp(argv[i+1], "3", 1) == 0 || strcmp(argv[i+1], "gff") == 0
              || strcmp(argv[i+1], "GFF") == 0)
        output = 3;
      else usage("Invalid output format specified.");
      i++;
    }
    else usage("Unknown option.");
  }

  /* Create the library context */
  ctx = prodigal_create(&config);
  if(ctx == NULL) {
    fprintf(stderr, "\nError: Failed to create Prodigal context.\n\n");
    exit(1);
  }

  /* Print header */
  if(quiet == 0) {
    fprintf(stderr, "-------------------------------------\n");
    fprintf(stderr, "PRODIGAL v%s [%s]         \n", VERSION, DATE);
    fprintf(stderr, "Univ of Tenn / Oak Ridge National Lab\n");
    fprintf(stderr, "Doug Hyatt, Loren Hauser, et al.     \n");
    fprintf(stderr, "-------------------------------------\n");
  }

  /* Read in the training file (if specified) */
  if(train_file != NULL) {
    if(config.meta_mode == 1) {
      fprintf(stderr, "\nError: cannot specify metagenomic sequence with a");
      fprintf(stderr, " training file.\n");
      exit(2);
    }
    rv = read_training_file(train_file, &ctx->tinf);
    if(rv == 1) do_training = 1;
    else {
      if(config.force_nonsd == 1) {
        fprintf(stderr, "\nError: cannot force non-SD finder with a training");
        fprintf(stderr, " file already created!\n"); exit(3);
      }
      if(quiet == 0)
        fprintf(stderr, "Reading in training data from file %s...", train_file);
      if(user_tt > 0 && user_tt != ctx->tinf.trans_table) {
        fprintf(stderr, "\n\nWarning: user-specified translation table does");
        fprintf(stderr, "not match the one in the specified training file! \n\n");
      }
      if(rv == -1) {
        fprintf(stderr, "\n\nError: training file did not read correctly!\n");
        exit(4);
      }
      ctx->trained = 1;
      if(quiet == 0) {
        fprintf(stderr, "done!\n");
        fprintf(stderr, "-------------------------------------\n");
      }
    }
  }

  /* Determine where standard input is coming from and react accordingly */
  if(config.meta_mode == 0 && train_file == NULL && input_file == NULL) {
    fnum = fileno(stdin);
    if(fstat(fnum, &fbuf) == -1) {
      fprintf(stderr, "\nError: can't fstat standard input.\n\n");
      exit(5);
    }
    if(S_ISCHR(fbuf.st_mode)) help();
    else if(S_ISREG(fbuf.st_mode)) { /* do nothing */ }
    else if(S_ISFIFO(fbuf.st_mode)) {
      piped = 1;
      if(copy_standard_input_to_file(input_copy, quiet) == -1) {
        fprintf(stderr, "\nError: can't copy stdin to file.\n\n");
        exit(5);
      }
      input_file = input_copy;
    }
  }

  /* Check i/o files (if specified) and prepare them for reading/writing */
  if(input_file != NULL) {
    input_ptr = INPUT_OPEN(input_file, "r");
    if(input_ptr == NULL) {
      fprintf(stderr, "\nError: can't open input file %s.\n\n", input_file);
      exit(5);
    }
  }
  if(input_ptr == NULL) {
    input_ptr = INPUT_OPEN("/dev/stdin", "r");
    if(input_ptr == NULL) {
      fprintf(stderr, "\nError: can't open input file %s.\n\n", input_file);
      exit(5);
    }
  }
  if(output_file != NULL) {
    output_ptr = fopen(output_file, "w");
    if(output_ptr == NULL) {
      fprintf(stderr, "\nError: can't open output file %s.\n\n", output_file);
      exit(6);
    }
  }
  if(start_file != NULL) {
    start_ptr = fopen(start_file, "w");
    if(start_ptr == NULL) {
      fprintf(stderr, "\nError: can't open start file %s.\n\n", start_file);
      exit(7);
    }
  }
  if(trans_file != NULL) {
    trans_ptr = fopen(trans_file, "w");
    if(trans_ptr == NULL) {
      fprintf(stderr, "\nError: can't open translation file %s.\n\n",
              trans_file);
      exit(8);
    }
  }
  if(nuc_file != NULL) {
    nuc_ptr = fopen(nuc_file, "w");
    if(nuc_ptr == NULL) {
      fprintf(stderr, "\nError: can't open gene nucleotide file %s.\n\n",
              nuc_file);
      exit(16);
    }
  }

  /***************************************************************************
    Single Genome Training:  Read in the sequence(s) and perform the
    training on them.  Uses the library context for all state.
  ***************************************************************************/
  if(config.meta_mode == 0 && (do_training == 1 || (do_training == 0 &&
     train_file == NULL))) {
    if(quiet == 0) {
      fprintf(stderr, "Request:  Single Genome, Phase:  Training\n");
      fprintf(stderr, "Reading in the sequence(s) to train...");
    }

    /* Read sequences directly into context buffers using existing FILE* I/O */
    ctx->slen = read_seq_training(input_ptr, ctx->seq, ctx->useq,
                                  &(ctx->tinf.gc), config.mask_regions,
                                  ctx->mlist, &ctx->nmask);
    if(ctx->slen == 0) {
      fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(9);
    }
    if(ctx->slen < MIN_SINGLE_GENOME) {
      fprintf(stderr, "\n\nError:  Sequence must be %d", MIN_SINGLE_GENOME);
      fprintf(stderr, " characters (only %d read).\n(Consider", ctx->slen);
      fprintf(stderr, " running with the -p meta option or finding");
      fprintf(stderr, " more contigs from the same genome.)\n\n");
      exit(10);
    }
    if(ctx->slen < IDEAL_SINGLE_GENOME) {
      fprintf(stderr, "\n\nWarning:  ideally Prodigal should be given at");
      fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
      fprintf(stderr, "training.\nYou may get better results with the ");
      fprintf(stderr, "-p meta option.\n\n");
    }
    rcom_seq(ctx->seq, ctx->rseq, ctx->useq, ctx->slen);
    ctx->gc = ctx->tinf.gc;
    if(quiet == 0) {
      fprintf(stderr, "%d bp seq created, %.2f pct GC\n", ctx->slen,
              ctx->tinf.gc*100.0);
    }

    /* Use library training pipeline */
    if(quiet == 0) {
      fprintf(stderr, "Locating all potential starts and stops...");
    }
    rv = prodigal_train(ctx);
    if(rv != PRODIGAL_OK) {
      fprintf(stderr, "\nError: training failed: %s\n",
              prodigal_last_error(ctx));
      exit(11);
    }
    if(quiet == 0) {
      fprintf(stderr, "done!\n");
    }

    /* If training specified, write the training file and exit. */
    if(do_training == 1) {
      if(quiet == 0) {
        fprintf(stderr, "Writing data to training file %s...", train_file);
      }
      rv = write_training_file(train_file, &ctx->tinf);
      if(rv != 0) {
        fprintf(stderr, "\nError: could not write training file!\n");
        exit(12);
      }
      else {
        if(quiet == 0) fprintf(stderr, "done!\n");
        prodigal_destroy(ctx);
        exit(0);
      }
    }

    /* Rewind input file */
    if(quiet == 0) fprintf(stderr, "-------------------------------------\n");
    if(INPUT_SEEK(input_ptr, 0, SEEK_SET) == -1) {
      fprintf(stderr, "\nError: could not rewind input file.\n");
      exit(13);
    }

    /* Reset sequence/dynamic programming variables */
    memset(ctx->seq, 0, (ctx->slen/4+1)*sizeof(unsigned char));
    memset(ctx->rseq, 0, (ctx->slen/4+1)*sizeof(unsigned char));
    memset(ctx->useq, 0, (ctx->slen/8+1)*sizeof(unsigned char));
    memset(ctx->nodes, 0, ctx->nn*sizeof(struct _node));
    ctx->nn = 0; ctx->slen = 0; ipath = 0; ctx->nmask = 0;
  }

  /* Initialize the training files for a metagenomic request */
  else if(config.meta_mode == 1) {
    if(quiet == 0) {
      fprintf(stderr, "Request:  Metagenomic, Phase:  Training\n");
      fprintf(stderr, "Initializing training files...");
    }
    /* Allocate and initialize metagenomic bins */
    ctx->meta = (struct _metagenomic_bin *)malloc(
        NUM_META * sizeof(struct _metagenomic_bin));
    if(ctx->meta == NULL) {
      fprintf(stderr, "\nError: Malloc failed on metagenomic bins.\n\n");
      exit(1);
    }
    for(i = 0; i < NUM_META; i++) {
      memset(&ctx->meta[i], 0, sizeof(struct _metagenomic_bin));
      strcpy(ctx->meta[i].desc, "None");
      ctx->meta[i].tinf = (struct _training *)malloc(sizeof(struct _training));
      if(ctx->meta[i].tinf == NULL) {
        fprintf(stderr, "\nError: Malloc failed on training structure.\n\n");
        exit(1);
      }
      memset(ctx->meta[i].tinf, 0, sizeof(struct _training));
    }
    initialize_metagenomic_bins(ctx->meta);
    ctx->meta_initialized = 1;
    if(quiet == 0) {
      fprintf(stderr, "done!\n");
      fprintf(stderr, "-------------------------------------\n");
    }
  }

  /* Print out header for gene finding phase */
  if(quiet == 0) {
    if(config.meta_mode == 1)
      fprintf(stderr, "Request:  Metagenomic, Phase:  Gene Finding\n");
    else fprintf(stderr, "Request:  Single Genome, Phase:  Gene Finding\n");
  }

  /* Read and process each sequence in the file in succession */
  sprintf(cur_header, "Prodigal_Seq_1");
  sprintf(new_header, "Prodigal_Seq_2");
  while((ctx->slen = next_seq_multi(input_ptr, ctx->seq, ctx->useq, &num_seq,
         &gc, config.mask_regions, ctx->mlist, &ctx->nmask, cur_header,
         new_header)) != -1) {
    rcom_seq(ctx->seq, ctx->rseq, ctx->useq, ctx->slen);
    if(ctx->slen == 0) {
      fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(14);
    }

    if(quiet == 0) {
      fprintf(stderr, "Finding genes in sequence #%d (%d bp)...", num_seq,
              ctx->slen);
    }

    /* Reallocate memory if this is the biggest sequence we've seen */
    if(ctx->slen > ctx->max_slen && ctx->slen > STT_NOD*8) {
      ctx->nodes = (struct _node *)realloc(ctx->nodes,
                   (int)(ctx->slen/8)*sizeof(struct _node));
      if(ctx->nodes == NULL) {
        fprintf(stderr, "Realloc failed on nodes\n\n");
        exit(11);
      }
      ctx->max_slen = ctx->slen;
    }

    /* Calculate short header for this sequence */
    calc_short_header(cur_header, short_header, num_seq);

    if(config.meta_mode == 0) { /* Single Genome Version */

      nn = add_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                     config.closed_ends, ctx->mlist, ctx->nmask, &ctx->tinf);
      qsort(ctx->nodes, nn, sizeof(struct _node), &compare_nodes);

      score_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, nn, &ctx->tinf,
                  config.closed_ends, config.meta_mode);
      if(start_ptr != stdout)
        write_start_file(start_ptr, ctx->nodes, nn, &ctx->tinf, num_seq,
                         ctx->slen, 0, NULL, VERSION, cur_header);
      record_overlapping_starts(ctx->nodes, nn, &ctx->tinf, 1);
      ipath = dprog(ctx->nodes, nn, &ctx->tinf, 1);
      eliminate_bad_genes(ctx->nodes, ipath, &ctx->tinf);
      ng = add_genes(ctx->genes, ctx->nodes, ipath);
      tweak_final_starts(ctx->genes, ng, ctx->nodes, nn, &ctx->tinf);
      record_gene_data(ctx->genes, ng, ctx->nodes, &ctx->tinf, num_seq);
      if(quiet == 0) {
        fprintf(stderr, "done!\n");
      }

      /* Output the genes */
      print_genes(output_ptr, ctx->genes, ng, ctx->nodes, ctx->slen, output,
                  num_seq, 0, NULL, &ctx->tinf, cur_header, short_header,
                  VERSION);
      fflush(output_ptr);
      if(trans_ptr != stdout)
        write_translations(trans_ptr, ctx->genes, ng, ctx->nodes, ctx->seq,
                           ctx->rseq, ctx->useq, ctx->slen, &ctx->tinf,
                           num_seq, short_header);
      if(nuc_ptr != stdout)
        write_nucleotide_seqs(nuc_ptr, ctx->genes, ng, ctx->nodes, ctx->seq,
                              ctx->rseq, ctx->useq, ctx->slen, &ctx->tinf,
                              num_seq, short_header);
    }

    else { /* Metagenomic Version */

      low = 0.88495*gc - 0.0102337;
      if(low > 0.65) low = 0.65;
      high = 0.86596*gc + .1131991;
      if(high < 0.35) high = 0.35;

      max_score = -100.0;
      for(i = 0; i < NUM_META; i++) {
        if(i == 0 || ctx->meta[i].tinf->trans_table !=
           ctx->meta[i-1].tinf->trans_table) {
          memset(ctx->nodes, 0, nn*sizeof(struct _node));
          nn = add_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                         config.closed_ends, ctx->mlist, ctx->nmask,
                         ctx->meta[i].tinf);
          qsort(ctx->nodes, nn, sizeof(struct _node), &compare_nodes);
        }
        if(ctx->meta[i].tinf->gc < low || ctx->meta[i].tinf->gc > high)
          continue;
        reset_node_scores(ctx->nodes, nn);
        score_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, nn,
                    ctx->meta[i].tinf, config.closed_ends, config.meta_mode);
        record_overlapping_starts(ctx->nodes, nn, ctx->meta[i].tinf, 1);
        ipath = dprog(ctx->nodes, nn, ctx->meta[i].tinf, 1);
        /* dprog returns -1 when it finds no gene path; upstream then reads
           nodes[-1].score, which is undefined (it crashes on macOS). The
           empty path scores 0.0 in dprog, where every node starts at 0, and
           that is also what upstream's read returns on glibc; see
           docs/upstream-behavior.md. */
        path_score = (ipath >= 0) ? ctx->nodes[ipath].score : 0.0;
        if(path_score > max_score) {
          max_phase = i;
          max_score = path_score;
          eliminate_bad_genes(ctx->nodes, ipath, ctx->meta[i].tinf);
          ng = add_genes(ctx->genes, ctx->nodes, ipath);
          tweak_final_starts(ctx->genes, ng, ctx->nodes, nn,
                             ctx->meta[i].tinf);
          record_gene_data(ctx->genes, ng, ctx->nodes, ctx->meta[i].tinf,
                           num_seq);
        }
      }

      /* Recover the nodes for the best of the runs */
      memset(ctx->nodes, 0, nn*sizeof(struct _node));
      nn = add_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                     config.closed_ends, ctx->mlist, ctx->nmask,
                     ctx->meta[max_phase].tinf);
      qsort(ctx->nodes, nn, sizeof(struct _node), &compare_nodes);
      score_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, nn,
                  ctx->meta[max_phase].tinf, config.closed_ends,
                  config.meta_mode);
      if(start_ptr != stdout)
        write_start_file(start_ptr, ctx->nodes, nn, ctx->meta[max_phase].tinf,
                         num_seq, ctx->slen, 1, ctx->meta[max_phase].desc,
                         VERSION, cur_header);

      if(quiet == 0) {
        fprintf(stderr, "done!\n");
      }

      /* Output the genes */
      print_genes(output_ptr, ctx->genes, ng, ctx->nodes, ctx->slen, output,
                  num_seq, 1, ctx->meta[max_phase].desc,
                  ctx->meta[max_phase].tinf, cur_header, short_header, VERSION);
      fflush(output_ptr);
      if(trans_ptr != stdout)
        write_translations(trans_ptr, ctx->genes, ng, ctx->nodes, ctx->seq,
                           ctx->rseq, ctx->useq, ctx->slen,
                           ctx->meta[max_phase].tinf, num_seq, short_header);
      if(nuc_ptr != stdout)
        write_nucleotide_seqs(nuc_ptr, ctx->genes, ng, ctx->nodes, ctx->seq,
                              ctx->rseq, ctx->useq, ctx->slen,
                              ctx->meta[max_phase].tinf, num_seq, short_header);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(ctx->seq, 0, (ctx->slen/4+1)*sizeof(unsigned char));
    memset(ctx->rseq, 0, (ctx->slen/4+1)*sizeof(unsigned char));
    memset(ctx->useq, 0, (ctx->slen/8+1)*sizeof(unsigned char));
    memset(ctx->nodes, 0, nn*sizeof(struct _node));
    nn = 0; ctx->slen = 0; ipath = 0; ctx->nmask = 0;
    strcpy(cur_header, new_header);
    sprintf(new_header, "Prodigal_Seq_%d\n", num_seq+1);
  }

  if(num_seq == 0) {
    fprintf(stderr, "\nError:  no input sequences to analyze.\n\n");
    exit(18);
  }

  /* Free all memory via library context */
  prodigal_destroy(ctx);

  /* Close all the filehandles and exit */
  INPUT_CLOSE(input_ptr);
  if(output_ptr != stdout) fclose(output_ptr);
  if(start_ptr != stdout) fclose(start_ptr);
  if(trans_ptr != stdout) fclose(trans_ptr);

  /* Remove tmp file */
  if(piped == 1 && remove(input_copy) != 0) {
    fprintf(stderr, "Could not delete tmp file %s.\n", input_copy);
    exit(18);
  }

  exit(0);
}

void version() {
  fprintf(stderr, "\nProdigal V%s: %s\n\n", VERSION, DATE);
  exit(0);
}

void usage(char *msg) {
  fprintf(stderr, "\n%s\n", msg);
  fprintf(stderr, "\nUsage:  prodigal [-a trans_file] [-c] [-d nuc_file]");
  fprintf(stderr, " [-f output_type]\n");
  fprintf(stderr, "                 [-g tr_table] [-h] [-i input_file] [-m]");
  fprintf(stderr, " [-n] [-o output_file]\n");
  fprintf(stderr, "                 [-p mode] [-q] [-s start_file]");
  fprintf(stderr, " [-t training_file] [-v]\n");
  fprintf(stderr, "\nDo 'prodigal -h' for more information.\n\n");
  exit(15);
}

void help() {
  fprintf(stderr, "\nUsage:  prodigal [-a trans_file] [-c] [-d nuc_file]");
  fprintf(stderr, " [-f output_type]\n");
  fprintf(stderr, "                 [-g tr_table] [-h] [-i input_file] [-m]");
  fprintf(stderr, " [-n] [-o output_file]\n");
  fprintf(stderr, "                 [-p mode] [-q] [-s start_file]");
  fprintf(stderr, " [-t training_file] [-v]\n");
  fprintf(stderr, "\n         -a:  Write protein translations to the selected ");
  fprintf(stderr, "file.\n");
  fprintf(stderr, "         -c:  Closed ends.  Do not allow genes to run off ");
  fprintf(stderr, "edges.\n");
  fprintf(stderr, "         -d:  Write nucleotide sequences of genes to the ");
  fprintf(stderr, "selected file.\n");
  fprintf(stderr, "         -f:  Select output format (gbk, gff, or sco).  ");
  fprintf(stderr, "Default is gbk.\n");
  fprintf(stderr, "         -g:  Specify a translation table to use (default");
  fprintf(stderr, " 11).\n");
  fprintf(stderr, "         -h:  Print help menu and exit.\n");
  fprintf(stderr, "         -i:  Specify FASTA/Genbank input file (default ");
  fprintf(stderr, "reads from stdin).\n");
  fprintf(stderr, "         -m:  Treat runs of N as masked sequence; don't");
  fprintf(stderr, " build genes across them.\n");
  fprintf(stderr, "         -n:  Bypass Shine-Dalgarno trainer and force");
  fprintf(stderr, " a full motif scan.\n");
  fprintf(stderr, "         -o:  Specify output file (default writes to ");
  fprintf(stderr, "stdout).\n");
  fprintf(stderr, "         -p:  Select procedure (single or meta).  Default");
  fprintf(stderr, " is single.\n");
  fprintf(stderr, "         -q:  Run quietly (suppress normal stderr output).\n");
  fprintf(stderr, "         -s:  Write all potential genes (with scores) to");
  fprintf(stderr, " the selected file.\n");
  fprintf(stderr, "         -t:  Write a training file (if none exists); ");
  fprintf(stderr, "otherwise, read and use\n");
  fprintf(stderr, "              the specified training file.\n");
  fprintf(stderr, "         -v:  Print version number and exit.\n\n");
  exit(0);
}

/* For piped input, we make a copy of stdin so we can rewind the file. */

int copy_standard_input_to_file(char *path, int quiet) {
  char line[MAX_LINE+1];
  FILE *wp;

  if(quiet == 0) {
    fprintf(stderr, "Piped input detected, copying stdin to a tmp file...");
  }

  wp = fopen(path, "w");
  if(wp == NULL) return -1;
  while(fgets(line, MAX_LINE, stdin) != NULL) {
    fprintf(wp, "%s", line);
  }
  fclose(wp);

  if(quiet == 0) {
    fprintf(stderr, "done!\n");
    fprintf(stderr, "-------------------------------------\n");
  }
  return 0;
}

#endif /* PRODIGAL_NO_MAIN */
