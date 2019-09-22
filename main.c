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
#include "sequence.h"
#include "metagenomic.h"
#include "node.h"
#include "dprog.h"
#include "gene.h"
#include "fptr.h"


#define VERSION "2.6.3"
#define DATE "February, 2016"

#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000


void version();
void usage(char *);
void help();
int copy_standard_input_to_file(char *, int);

int main(int argc, char *argv[]) {

  int rv, slen, nn, ng, i, ipath, *gc_frame, do_training, output, max_phase;
  int closed, do_mask, nmask, force_nonsd, user_tt, is_meta, num_seq, quiet;
  int piped, max_slen, fnum;
  double max_score, gc, low, high;
  unsigned char *seq, *rseq, *useq;
  char *train_file, *start_file, *trans_file, *nuc_file; 
  char *input_file, *output_file, input_copy[MAX_LINE];
  char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
  FILE *output_ptr, *start_ptr, *trans_ptr, *nuc_ptr;
  fptr input_ptr = NULL;
  struct stat fbuf;
  pid_t pid;
  struct _node *nodes;
  struct _gene *genes;
  struct _training tinf;
  struct _metagenomic_bin meta[NUM_META];
  mask mlist[MAX_MASKS];

  /* Allocate memory and initialize variables */
  seq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char));
  rseq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char));
  useq = (unsigned char *)malloc(MAX_SEQ/8*sizeof(unsigned char));
  nodes = (struct _node *)malloc(STT_NOD*sizeof(struct _node));
  genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene));
  if(seq == NULL || rseq == NULL || nodes == NULL || genes == NULL) {
    fprintf(stderr, "\nError: Malloc failed on sequence/orfs\n\n"); exit(1);
  }
  memset(seq, 0, MAX_SEQ/4*sizeof(unsigned char));
  memset(rseq, 0, MAX_SEQ/4*sizeof(unsigned char));
  memset(useq, 0, MAX_SEQ/8*sizeof(unsigned char));
  memset(nodes, 0, STT_NOD*sizeof(struct _node));
  memset(genes, 0, MAX_GENES*sizeof(struct _gene));
  memset(&tinf, 0, sizeof(struct _training));

  for(i = 0; i < NUM_META; i++) {
    memset(&meta[i], 0, sizeof(struct _metagenomic_bin));
    strcpy(meta[i].desc, "None");
    meta[i].tinf = (struct _training *)malloc(sizeof(struct _training));
    if(meta[i].tinf == NULL) {
      fprintf(stderr, "\nError: Malloc failed on training structure.\n\n"); 
      exit(1);
    }
    memset(meta[i].tinf, 0, sizeof(struct _training));
  }
  nn = 0; slen = 0; ipath = 0; ng = 0; nmask = 0;
  user_tt = 0; is_meta = 0; num_seq = 0; quiet = 0;
  max_phase = 0; max_score = -100.0;
  train_file = NULL; do_training = 0;
  start_file = NULL; trans_file = NULL; nuc_file = NULL;
  start_ptr = stdout; trans_ptr = stdout; nuc_ptr = stdout;
  input_file = NULL; output_file = NULL; piped = 0;
  output_ptr = stdout; max_slen = 0;
  output = 0; closed = 0; do_mask = 0; force_nonsd = 0;

  /* Filename for input copy if needed */
  pid = getpid();
  sprintf(input_copy, "tmp.prodigal.stdin.%d", pid);

  /***************************************************************************
    Set the start score weight.  Changing this number can dramatically
    affect the performance of the program.  Some genomes want it high (6+),
    and some prefer it low (2.5-3).  Attempts were made to determine this 
    weight dynamically, but none were successful.  Therefore, we just 
    manually set the weight to an average value that seems to work decently 
    for 99% of genomes.  This problem may be revisited in future versions.
  ***************************************************************************/
  tinf.st_wt = 4.35;
  tinf.trans_table = 11;

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
      closed = 1;
    else if(strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "-Q") == 0)
      quiet = 1;
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-M") == 0)
      do_mask = 1;
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "-N") == 0)
      force_nonsd = 1;
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
      tinf.trans_table = atoi(argv[i+1]);
      if(tinf.trans_table < 1 || tinf.trans_table > 25 || tinf.trans_table == 7
         || tinf.trans_table == 8 || (tinf.trans_table >= 17 && tinf.trans_table
         <= 20))
        usage("Invalid translation table specified.");
      user_tt = tinf.trans_table;
      i++;
    }
    else if(strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-P") == 0) {
      if(argv[i+1][0] == '0' || argv[i+1][0] == 's' || argv[i+1][0] ==
              'S') is_meta = 0;
      else if(argv[i+1][0] == '1' || argv[i+1][0] == 'm' || argv[i+1][0] ==
              'M') is_meta = 1; 
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
    if(is_meta == 1) {
      fprintf(stderr, "\nError: cannot specify metagenomic sequence with a");
      fprintf(stderr, " training file.\n");
      exit(2);
    } 
    rv = read_training_file(train_file, &tinf);
    if(rv == 1) do_training = 1;
    else {
      if(force_nonsd == 1) { 
        fprintf(stderr, "\nError: cannot force non-SD finder with a training");
        fprintf(stderr, " file already created!\n"); exit(3);
      }
      if(quiet == 0)
        fprintf(stderr, "Reading in training data from file %s...", train_file);
      if(user_tt > 0 && user_tt != tinf.trans_table) { 
        fprintf(stderr, "\n\nWarning: user-specified translation table does");
        fprintf(stderr, "not match the one in the specified training file! \n\n");
      }
      if(rv == -1) { 
        fprintf(stderr, "\n\nError: training file did not read correctly!\n"); 
        exit(4); 
      }
      if(quiet == 0) {
        fprintf(stderr, "done!\n"); 
        fprintf(stderr, "-------------------------------------\n");
      }
    }
  }

  /* Determine where standard input is coming from and react accordingly */
  if(is_meta == 0 && train_file == NULL && input_file == NULL) {
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
    training on them.
  ***************************************************************************/
  if(is_meta == 0 && (do_training == 1 || (do_training == 0 && train_file == 
     NULL))) {
    if(quiet == 0) {
      fprintf(stderr, "Request:  Single Genome, Phase:  Training\n");
      fprintf(stderr, "Reading in the sequence(s) to train..."); 
    }
    slen = read_seq_training(input_ptr, seq, useq, &(tinf.gc), do_mask, mlist,
                             &nmask);
    if(slen == 0) {
      fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(9);
    }
    if(slen < MIN_SINGLE_GENOME) {
      fprintf(stderr, "\n\nError:  Sequence must be %d", MIN_SINGLE_GENOME);
      fprintf(stderr, " characters (only %d read).\n(Consider", slen);
      fprintf(stderr, " running with the -p meta option or finding");
      fprintf(stderr, " more contigs from the same genome.)\n\n");
      exit(10);
    }
    if(slen < IDEAL_SINGLE_GENOME) {
      fprintf(stderr, "\n\nWarning:  ideally Prodigal should be given at");
      fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
      fprintf(stderr, "training.\nYou may get better results with the ");
      fprintf(stderr, "-p meta option.\n\n");
    }
    rcom_seq(seq, rseq, useq, slen);
    if(quiet == 0) {
      fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc*100.0);
    }

    /***********************************************************************
      Find all the potential starts and stops, sort them, and create a 
      comprehensive list of nodes for dynamic programming.
    ***********************************************************************/
    if(quiet == 0) {
      fprintf(stderr, "Locating all potential starts and stops..."); 
    }
    if(slen > max_slen && slen > STT_NOD*8) {
      nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
      if(nodes == NULL) {
        fprintf(stderr, "Realloc failed on nodes\n\n");
        exit(11);
      }
      max_slen = slen;
    }
    nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
    qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
    if(quiet == 0) {
      fprintf(stderr, "%d nodes\n", nn); 
    }

    /***********************************************************************
      Scan all the ORFS looking for a potential GC bias in a particular
      codon position.  This information will be used to acquire a good
      initial set of genes.
    ***********************************************************************/
    if(quiet == 0) {
      fprintf(stderr, "Looking for GC bias in different frames...");
    }
    gc_frame = calc_most_gc_frame(seq, slen);
    if(gc_frame == NULL) {
      fprintf(stderr, "Malloc failed on gc frame plot\n\n");
      exit(11);
    }
    record_gc_bias(gc_frame, nodes, nn, &tinf);
    if(quiet == 0) {
      fprintf(stderr, "frame bias scores: %.2f %.2f %.2f\n", tinf.bias[0],
              tinf.bias[1], tinf.bias[2]); 
    }
    free(gc_frame);

    /***********************************************************************
      Do an initial dynamic programming routine with just the GC frame
      bias used as a scoring function.  This will get an initial set of 
      genes to train on. 
    ***********************************************************************/
    if(quiet == 0) {
      fprintf(stderr, "Building initial set of genes to train from...");
    }
    record_overlapping_starts(nodes, nn, &tinf, 0);
    ipath = dprog(nodes, nn, &tinf, 0);
    if(quiet == 0) {
      fprintf(stderr, "done!\n"); 
    }

    /***********************************************************************
      Gather dicodon statistics for the training set.  Score the entire set
      of nodes.                               
    ***********************************************************************/
    if(quiet == 0) {
      fprintf(stderr, "Creating coding model and scoring nodes...");
    }
    calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);
    if(quiet == 0) {
      fprintf(stderr, "done!\n"); 
    }

    /***********************************************************************
      Determine if this organism uses Shine-Dalgarno or not and score the 
      nodes appropriately.
    ***********************************************************************/
    if(quiet == 0) {
      fprintf(stderr, "Examining upstream regions and training starts...");
    }
    rbs_score(seq, rseq, slen, nodes, nn, &tinf);
    train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
    determine_sd_usage(&tinf);
    if(force_nonsd == 1) tinf.uses_sd = 0;
    if(tinf.uses_sd == 0) train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf);
    if(quiet == 0) {
      fprintf(stderr, "done!\n"); 
    }

    /* If training specified, write the training file and exit. */
    if(do_training == 1) {
      if(quiet == 0) {
        fprintf(stderr, "Writing data to training file %s...", train_file);
      }
      rv = write_training_file(train_file, &tinf);
      if(rv != 0) { 
        fprintf(stderr, "\nError: could not write training file!\n"); 
        exit(12); 
      }
      else { 
        if(quiet == 0) fprintf(stderr, "done!\n"); 
        exit(0); 
      }
    }

    /* Rewind input file */    
    if(quiet == 0) fprintf(stderr, "-------------------------------------\n");
    if(INPUT_SEEK(input_ptr, 0, SEEK_SET) == -1) {
      fprintf(stderr, "\nError: could not rewind input file.\n"); 
      exit(13);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;
  }

  /* Initialize the training files for a metagenomic request */
  else if(is_meta == 1) {
    if(quiet == 0) {
      fprintf(stderr, "Request:  Metagenomic, Phase:  Training\n");
      fprintf(stderr, "Initializing training files...");
    }
    initialize_metagenomic_bins(meta);
    if(quiet == 0) {
      fprintf(stderr, "done!\n");
      fprintf(stderr, "-------------------------------------\n");
    }
  }

  /* Print out header for gene finding phase */
  if(quiet == 0) {
    if(is_meta == 1) 
      fprintf(stderr, "Request:  Metagenomic, Phase:  Gene Finding\n");
    else fprintf(stderr, "Request:  Single Genome, Phase:  Gene Finding\n");
  }

  /* Read and process each sequence in the file in succession */
  sprintf(cur_header, "Prodigal_Seq_1");
  sprintf(new_header, "Prodigal_Seq_2");
  while((slen = next_seq_multi(input_ptr, seq, useq, &num_seq, &gc, 
         do_mask, mlist, &nmask, cur_header, new_header)) != -1) {
    rcom_seq(seq, rseq, useq, slen);
    if(slen == 0) {
      fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(14);
    }

    if(quiet == 0) {
      fprintf(stderr, "Finding genes in sequence #%d (%d bp)...", num_seq, slen);
    }

    /* Reallocate memory if this is the biggest sequence we've seen */
    if(slen > max_slen && slen > STT_NOD*8) {
      nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
      if(nodes == NULL) {
        fprintf(stderr, "Realloc failed on nodes\n\n");
        exit(11);
      }
      max_slen = slen;
    }

    /* Calculate short header for this sequence */
    calc_short_header(cur_header, short_header, num_seq);

    if(is_meta == 0) { /* Single Genome Version */

      /***********************************************************************
        Find all the potential starts and stops, sort them, and create a 
        comprehensive list of nodes for dynamic programming.
      ***********************************************************************/
      nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

      /***********************************************************************
        Second dynamic programming, using the dicodon statistics as the
        scoring function.                                
      ***********************************************************************/
      score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, is_meta);
      if(start_ptr != stdout) 
        write_start_file(start_ptr, nodes, nn, &tinf, num_seq, slen, 0, NULL,
                         VERSION, cur_header);
      record_overlapping_starts(nodes, nn, &tinf, 1);
      ipath = dprog(nodes, nn, &tinf, 1);
      eliminate_bad_genes(nodes, ipath, &tinf);
      ng = add_genes(genes, nodes, ipath);
      tweak_final_starts(genes, ng, nodes, nn, &tinf);
      record_gene_data(genes, ng, nodes, &tinf, num_seq);
      if(quiet == 0) {
        fprintf(stderr, "done!\n"); 
      }

      /* Output the genes */
      print_genes(output_ptr, genes, ng, nodes, slen, output, num_seq, 0, NULL,
                  &tinf, cur_header, short_header, VERSION);
      fflush(output_ptr);
      if(trans_ptr != stdout)
        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                              &tinf, num_seq, short_header);
      if(nuc_ptr != stdout)
        write_nucleotide_seqs(nuc_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                              &tinf, num_seq, short_header);
    }

    else { /* Metagenomic Version */

      low = 0.88495*gc - 0.0102337;
      if(low > 0.65) low = 0.65;
      high = 0.86596*gc + .1131991;
      if(high < 0.35) high = 0.35;

      max_score = -100.0;
      for(i = 0; i < NUM_META; i++) { 
        if(i == 0 || meta[i].tinf->trans_table != 
           meta[i-1].tinf->trans_table) {
          memset(nodes, 0, nn*sizeof(struct _node));
          nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, 
                         meta[i].tinf);
          qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        }
        if(meta[i].tinf->gc < low || meta[i].tinf->gc > high) continue;  
        reset_node_scores(nodes, nn);
        score_nodes(seq, rseq, slen, nodes, nn, meta[i].tinf, closed, is_meta);
        record_overlapping_starts(nodes, nn, meta[i].tinf, 1);
        ipath = dprog(nodes, nn, meta[i].tinf, 1);
        if(nodes[ipath].score > max_score) {
          max_phase = i;
          max_score = nodes[ipath].score;
          eliminate_bad_genes(nodes, ipath, meta[i].tinf);
          ng = add_genes(genes, nodes, ipath);
          tweak_final_starts(genes, ng, nodes, nn, meta[i].tinf);
          record_gene_data(genes, ng, nodes, meta[i].tinf, num_seq);
        }
      }    

      /* Recover the nodes for the best of the runs */
      memset(nodes, 0, nn*sizeof(struct _node));
      nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                       meta[max_phase].tinf);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
      score_nodes(seq, rseq, slen, nodes, nn, meta[max_phase].tinf, closed,
                  is_meta);
      if(start_ptr != stdout) 
        write_start_file(start_ptr, nodes, nn, meta[max_phase].tinf, 
                         num_seq, slen, 1, meta[max_phase].desc, VERSION,
                         cur_header);

      if(quiet == 0) {
        fprintf(stderr, "done!\n"); 
      }

      /* Output the genes */
      print_genes(output_ptr, genes, ng, nodes, slen, output, num_seq, 1,
                  meta[max_phase].desc, meta[max_phase].tinf, cur_header, 
                  short_header, VERSION);
      fflush(output_ptr);
      if(trans_ptr != stdout)
        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                           meta[max_phase].tinf, num_seq, short_header);
      if(nuc_ptr != stdout)
        write_nucleotide_seqs(nuc_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                              meta[max_phase].tinf, num_seq, short_header);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;
    strcpy(cur_header, new_header);
    sprintf(new_header, "Prodigal_Seq_%d\n", num_seq+1);
  }

  if(num_seq == 0) {
    fprintf(stderr, "\nError:  no input sequences to analyze.\n\n");
    exit(18);
  }

  /* Free all memory */
  free(seq);
  free(rseq);
  free(useq);
  free(nodes);
  free(genes);
  for(i = 0; i < NUM_META; i++) free(meta[i].tinf);

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
