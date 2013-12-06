/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2013 University of Tennessee / UT-Battelle

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
#include <ctype.h>
#include "sequence.h"
#include "anonymous.h"
#include "node.h"
#include "dprog.h"
#include "gene.h"

#define VERSION "3.0-devel/Nov2013"
#define DATE "December, 2013"

#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000

void version();
void usage(char *);
void help();
int copy_standard_input_to_file(char *, int);

int main(int argc, char *argv[]) {

  int rv, slen, nn, ng, i, j, ipath, *gc_frame, outfmt, max_phase;
  int closed, num_seq, cross_gaps, quiet;
  int piped, max_slen, fnum;
  int mode; /* 0 = normal, 1 = training, 2 = anonymous/metagenomic */
  int force_nonsd = 0, is_anon = 0; /* deprecated and slated for removal PDH */
  double max_score, gc, low, high;
  unsigned char *seq, *rseq, *useq;
  char *train_file, *start_file, *trans_file, *nuc_file; 
  char *input_file, *output_file, input_copy[MAX_LINE];
  char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
  FILE *input_ptr, *output_ptr, *start_ptr, *trans_ptr, *nuc_ptr;
  struct stat fbuf;
  pid_t pid;
  struct _node *nodes;
  struct _gene *genes;
  struct _training tinf;
  struct _preset_genome_bin presets[NUM_PRESET_GENOME];

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

  for(i = 0; i < NUM_PRESET_GENOME; i++) {
    memset(&presets[i], 0, sizeof(struct _preset_genome_bin));
    strcpy(presets[i].desc, "None");
    presets[i].tinf = (struct _training *)malloc(sizeof(struct _training));
    if(presets[i].tinf == NULL) {
      fprintf(stderr, "\nError: Malloc failed on training structure.\n\n"); 
      exit(1);
    }
    memset(presets[i].tinf, 0, sizeof(struct _training));
  }
  nn = 0; slen = 0; ipath = 0; ng = 0; cross_gaps = 0;
  mode = 0; num_seq = 0; quiet = 0;
  max_phase = 0; max_score = -100.0;
  train_file = NULL;
  start_file = NULL; trans_file = NULL; nuc_file = NULL;
  start_ptr = stdout; trans_ptr = stdout; nuc_ptr = stdout;
  input_file = NULL; output_file = NULL; piped = 0;
  input_ptr = stdin; output_ptr = stdout; max_slen = 0;
  outfmt = -1; closed = 0;

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
  tinf.trans_table = -1;

  /* Parse the command line arguments */
  for(i = 1; i < argc; i++) {
    if(argv[i][0] == '-')
      for(j = 0; j < strlen(argv[i]); j++) 
        argv[i][j] = tolower(argv[i][j]);
  }
  for(i = 1; i < argc; i++) {
    if((i == argc-1 || argv[i+1][0] == '-') && 
       (strcmp(argv[i], "-t") == 0 || 
       strcmp(argv[i], "--training_file") == 0 ||
       strcmp(argv[i], "-a") == 0 || 
       strcmp(argv[i], "--protein_file") == 0 ||
       strcmp(argv[i], "-d") == 0 || 
       strcmp(argv[i], "--mrna_file") == 0 ||
       strcmp(argv[i], "-g") == 0 || 
       strcmp(argv[i], "--trans_table") == 0 ||
       strcmp(argv[i], "-f") == 0 || 
       strcmp(argv[i], "--output_format") == 0 ||
       strcmp(argv[i], "-s") == 0 || 
       strcmp(argv[i], "--start_file") == 0 ||
       strcmp(argv[i], "-i") == 0 || 
       strcmp(argv[i], "--input_file") == 0 ||
       strcmp(argv[i], "-o") == 0 || 
       strcmp(argv[i], "--output_file") == 0 ||
       strcmp(argv[i], "-m") == 0 || 
       strcmp(argv[i], "--mode") == 0))
      usage("-a/-d/-f/-g/-i/-m/-o/-s/-t options require valid parameters.");
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--nopartial") == 0)
      closed = 1;
    else if(strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--quiet") == 0)
      quiet = 1;
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nogaps") == 0)
      cross_gaps = 1;
    else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) 
      help();
    else if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) 
      version();
    else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--protein_file") 
            == 0) {
      trans_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--mrna_file") == 0) {
      nuc_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-i") == 0 || 
            strcmp(argv[i], "--input_file") == 0) {
      input_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || 
            strcmp(argv[i], "--output_file") == 0) {
      output_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || 
            strcmp(argv[i], "--start_file") == 0) {
      start_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || 
            strcmp(argv[i], "--training_file") == 0) {
      train_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-g") == 0 || 
            strcmp(argv[i], "--trans_table") == 0) {
      tinf.trans_table = atoi(argv[i+1]);
      if(tinf.trans_table < 0 || tinf.trans_table > 25 || tinf.trans_table == 7
         || tinf.trans_table == 8 || (tinf.trans_table >= 17 && tinf.trans_table
         <= 20))
        usage("Invalid translation table specified.");
      i++;
    }
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mode") == 0) {
      if(argv[i+1][0] == 'n') mode = 0;
      else if(argv[i+1][0] == 't') mode = 1; 
      else if(argv[i+1][0] == 'a') mode = 2; 
      else usage("Invalid mode specified (should be normal, train, or anon).");
      i++;
    }
    else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--output_format") 
            == 0) {
      if(strncmp(argv[i+1], "0", 1) == 0 || strcmp(argv[i+1], "gbk") == 0)
        outfmt = 0;
      else if(strncmp(argv[i+1], "1", 1) == 0 || strcmp(argv[i+1], "gca") == 0)
        outfmt = 1;
      else if(strncmp(argv[i+1], "2", 1) == 0 || strcmp(argv[i+1], "sco") == 0)
        outfmt = 2;
      else if(strncmp(argv[i+1], "3", 1) == 0 || strcmp(argv[i+1], "gff") == 0)
        outfmt = 3;
      else if(strncmp(argv[i+1], "4", 1) == 0 || strcmp(argv[i+1], "sqn") == 0)
        outfmt = 4;
      else usage("Invalid output format specified.");
      i++;
    }
    else {
      fprintf(stderr, "Unknown option '%s'.", argv[i]);
      usage("");
    }
  }

  /* Training mode can't have output format or extra files specified. */
  if(mode == 1 && (start_file != NULL || nuc_file != NULL ||
     trans_file != NULL || outfmt != -1)) {
    usage("-a/-d/-f/-s options cannot be used in training mode.");
  }

  /* Normal/anonymous can't have training files specified. */
  if(mode != 0 && train_file != NULL) {
    usage("Can only specify training file in normal mode.");
  }

  /* Anonymous mode can't have a specified value for genetic code. */
  /* Nor can normal mode if using a training file. */
  if((mode == 2 || (mode == 0 && train_file != NULL)) && 
     tinf.trans_table != -1) {
    usage("Can't specify translation table with anonymous mode or a training file.");
  }

  if(outfmt == -1) outfmt = 3; /* gff default output format */

  /* Print header */
  if(quiet == 0) {
    fprintf(stderr, "-------------------------------------\n");
    fprintf(stderr, "PRODIGAL v%s [%s]         \n", VERSION, DATE);
    fprintf(stderr, "Univ of Tenn / Oak Ridge National Lab\n");
    fprintf(stderr, "Doug Hyatt, Loren Hauser, et al.     \n");
    fprintf(stderr, "-------------------------------------\n");
  }
  if(quiet == 0) {
    if(mode == 0) fprintf(stderr, "Mode: Normal, Phase: Training\n");
    else if(mode == 1) fprintf(stderr, "Mode: Training, Phase: Training\n");
    else if(mode == 2) fprintf(stderr, "Mode: Anonymous, Phase: Training\n");
  }
  
  /* If we're in normal mode and not reading from a training file, */
  /* then we have to make two passes over the sequence.  Since we  */
  /* rewind after the first pass (something Windows can't do to    */
  /* stdin), we copy stdin to a temp file so that we can rewind it */
  /* in Windows. If there's nothing present on stdin, we print     */
  /* the help message and exit. */
  if(mode == 0 && train_file == NULL && input_file == NULL) {
    fnum = fileno(stdin);
    if(fstat(fnum, &fbuf) == -1) {
      fprintf(stderr, "\nError: can't fstat standard input.\n\n");
      exit(4);
    }
    if(S_ISCHR(fbuf.st_mode)) {
      if(argc == 1) help();
      else {
        fprintf(stderr, "\nError: options specified but no input ");
        fprintf(stderr, "detected.\n\n");
        exit(4);
      }
    }
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

  /* Read in the training file (if specified) */
  if(train_file != NULL) {
    if(quiet == 0)
      fprintf(stderr, "Reading in training data from file %s...", train_file);
    rv = read_training_file(train_file, &tinf);
    if(rv == -1) { 
      fprintf(stderr, "\n\nError: training file did not read correctly!\n"); 
      exit(2); 
    }
    if(quiet == 0) {
      fprintf(stderr, "done!\n"); 
      fprintf(stderr, "-------------------------------------\n");
    }
  }

  /* Check i/o files (if specified) and prepare them for reading/writing */
  if(input_file != NULL) {
    input_ptr = fopen(input_file, "r");
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
  if(mode == 1 || (mode == 0 && train_file == NULL)) {
    if(quiet == 0) {
      fprintf(stderr, "Reading in the sequence(s) to train..."); 
    }
    slen = read_seq_training(input_ptr, seq, useq, &(tinf.gc), closed);
    if(slen == 0) {
      fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
      fprintf(stderr, "Genbank, or EMBL format).\n\n");
      exit(9);
    }
    if(slen < MIN_SINGLE_GENOME) {
      fprintf(stderr, "\n\nError:  Sequence must be %d", MIN_SINGLE_GENOME);
      fprintf(stderr, " characters (only %d read).\n(Consider", slen);
      fprintf(stderr, " running with the -m anon option or finding");
      fprintf(stderr, " more contigs from the same genome.)\n\n");
      exit(10);
    }
    if(slen < IDEAL_SINGLE_GENOME) {
      fprintf(stderr, "\n\nWarning:  ideally Prodigal should be given at");
      fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
      fprintf(stderr, "training.\nYou may get better results with the ");
      fprintf(stderr, "'-m anon' option.\n\n");
    }
    rcom_seq(seq, rseq, useq, slen);
    if(quiet == 0) {
      fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc*100.0);
    }

    /***********************************************************************
      Scan the sequence and deduce what translation table it uses.
    ***********************************************************************/
    if(tinf.trans_table == -1) {
      if(quiet == 0) fprintf(stderr, "Detecting translation table...");
      tinf.trans_table = detect_translation_table(seq, rseq, useq, slen, tinf.gc);
      if(quiet == 0) fprintf(stderr, "table %d selected\n", tinf.trans_table);
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
    nn = add_nodes(seq, rseq, slen, nodes, closed, cross_gaps, &tinf);
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
    if(quiet == 0) fprintf(stderr, "done!\n");

    /***********************************************************************
      Gather dicodon statistics for the training set.  Score the entire set
      of nodes.                               
    ***********************************************************************/
    if(quiet == 0) {
      fprintf(stderr, "Creating coding model and scoring nodes...");
    }
    calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);
    if(quiet == 0) fprintf(stderr, "done!\n");

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
    if(quiet == 0) fprintf(stderr, "done!\n");

    /* If training specified, write the training file and exit. */
    if(mode == 1) {
      if(quiet == 0) {
        fprintf(stderr, "Writing data to training file %s...", output_file);
      }
      rv = write_training_file(output_file, &tinf);
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
    if(fseek(input_ptr, 0, SEEK_SET) == -1) {
      fprintf(stderr, "\nError: could not rewind input file.\n"); 
      exit(13);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0;
  }

  /* Initialize the training files for an anonymous request */
  else if(mode == 2) {
    if(quiet == 0) {
      fprintf(stderr, "Initializing preset training files...");
    }
    initialize_preset_genome_bins(presets);
    if(quiet == 0) {
      fprintf(stderr, "done!\n");
      fprintf(stderr, "-------------------------------------\n");
    }
  }

  /* Print out header for gene finding phase */
  if(quiet == 0) {
    if(mode == 2) 
      fprintf(stderr, "Mode: Anonymous, Phase: Gene Finding\n");
    else fprintf(stderr, "Mode: Normal, Phase: Gene Finding\n");
  }

  /* Read and process each sequence in the file in succession */
  sprintf(cur_header, "Prodigal_Seq_1");
  sprintf(new_header, "Prodigal_Seq_2");
  while((slen = next_seq_multi(input_ptr, seq, useq, &num_seq, &gc, 
         cur_header, new_header)) != -1) {
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

    if(mode != 2) { /* Single Genome Version */

      /***********************************************************************
        Find all the potential starts and stops, sort them, and create a 
        comprehensive list of nodes for dynamic programming.
      ***********************************************************************/
      nn = add_nodes(seq, rseq, slen, nodes, closed, cross_gaps, &tinf);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

      /***********************************************************************
        Second dynamic programming, using the dicodon statistics as the
        scoring function.                                
      ***********************************************************************/
      score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, is_anon);
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
      print_genes(output_ptr, genes, ng, nodes, slen, outfmt, num_seq, 0, NULL,
                  &tinf, cur_header, short_header, VERSION);
      fflush(output_ptr);
      if(trans_ptr != stdout)
        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                              &tinf, num_seq, short_header);
      if(nuc_ptr != stdout)
        write_nucleotide_seqs(nuc_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                              &tinf, num_seq, short_header);
    }

    else { /* Anonymous (Metagenomic) Version */
is_anon = 1;  /* deprecated slated for removal PDH */
      low = 0.88495*gc - 0.0102337;
      if(low > 0.65) low = 0.65;
      high = 0.86596*gc + .1131991;
      if(high < 0.35) high = 0.35;

      max_score = -100.0;
      for(i = 0; i < NUM_PRESET_GENOME; i++) { 
        if(i == 0 || presets[i].tinf->trans_table != 
           presets[i-1].tinf->trans_table) {
          memset(nodes, 0, nn*sizeof(struct _node));
          nn = add_nodes(seq, rseq, slen, nodes, closed, cross_gaps, presets[i].tinf);
          qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        }
        if(presets[i].tinf->gc < low || presets[i].tinf->gc > high) continue;  
        reset_node_scores(nodes, nn);
        score_nodes(seq, rseq, slen, nodes, nn, presets[i].tinf, closed, is_anon);
        record_overlapping_starts(nodes, nn, presets[i].tinf, 1);
        ipath = dprog(nodes, nn, presets[i].tinf, 1);
        if(nodes[ipath].score > max_score) {
          max_phase = i;
          max_score = nodes[ipath].score;
          eliminate_bad_genes(nodes, ipath, presets[i].tinf);
          ng = add_genes(genes, nodes, ipath);
          tweak_final_starts(genes, ng, nodes, nn, presets[i].tinf);
          record_gene_data(genes, ng, nodes, presets[i].tinf, num_seq);
        }
      }    

      /* Recover the nodes for the best of the runs */
      memset(nodes, 0, nn*sizeof(struct _node));
      nn = add_nodes(seq, rseq, slen, nodes, closed, cross_gaps, presets[max_phase].tinf);
      qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
      score_nodes(seq, rseq, slen, nodes, nn, presets[max_phase].tinf, closed,
                  is_anon);
      if(start_ptr != stdout) 
        write_start_file(start_ptr, nodes, nn, presets[max_phase].tinf, 
                         num_seq, slen, 1, presets[max_phase].desc, VERSION,
                         cur_header);

      if(quiet == 0) {
        fprintf(stderr, "done!\n"); 
      }

      /* Output the genes */
      print_genes(output_ptr, genes, ng, nodes, slen, outfmt, num_seq, 1,
                  presets[max_phase].desc, presets[max_phase].tinf, cur_header, 
                  short_header, VERSION);
      fflush(output_ptr);
      if(trans_ptr != stdout)
        write_translations(trans_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                           presets[max_phase].tinf, num_seq, short_header);
      if(nuc_ptr != stdout)
        write_nucleotide_seqs(nuc_ptr, genes, ng, nodes, seq, rseq, useq, slen,
                              presets[max_phase].tinf, num_seq, short_header);
    }

    /* Reset all the sequence/dynamic programming variables */
    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0;
    strcpy(cur_header, new_header);
    sprintf(new_header, "Prodigal_Seq_%d\n", num_seq+1);
  }

  if(num_seq == 0) {
    fprintf(stderr, "\nError:  no input sequences to analyze.\n\n");
    exit(18);
  }

  /* Free all memory */
  if(seq != NULL) free(seq);
  if(rseq != NULL) free(rseq);
  if(useq != NULL) free(useq);
  if(nodes != NULL) free(nodes);
  if(genes != NULL) free(genes);
  for(i = 0; i < NUM_PRESET_GENOME; i++) if(presets[i].tinf != NULL) free(presets[i].tinf);

  /* Close all the filehandles and exit */
  if(input_ptr != stdin) fclose(input_ptr);
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
  fprintf(stderr, "\nError: %s\n", msg);
  fprintf(stderr, "\nUsage:  prodigal [-a protein_file] [-c] [-d mrna_file]");
  fprintf(stderr, " [-f out_format]\n");
  fprintf(stderr, "                 [-g tr_table] [-h] [-i input_file]");
  fprintf(stderr, " [-m mode] [-n]\n");
  fprintf(stderr, "                 [-o output_file] [-q] [-s start_file]");
  fprintf(stderr, " [-t train_file]\n                 [-v]\n");
  fprintf(stderr, "\nDo 'prodigal -h' for more information.\n\n");
  exit(15);
}

void help() {
  fprintf(stderr, "\nUsage:  prodigal [-a protein_file] [-c] [-d mrna_file]");
  fprintf(stderr, " [-f out_format]\n");
  fprintf(stderr, "                 [-g tr_table] [-h] [-i input_file]");
  fprintf(stderr, " [-m mode] [-n]\n");
  fprintf(stderr, "                 [-o output_file] [-q] [-s start_file]");
  fprintf(stderr, " [-t train_file]\n                 [-v]\n");
  fprintf(stderr, "\nGene Modeling Parameters\n\n");
  fprintf(stderr, "  -m, --mode:           Specify mode (normal, train, or anon).\n");
  fprintf(stderr, "                          normal:   Single genome, any number of\n");
  fprintf(stderr, "                                    sequences. (Default)\n");
  fprintf(stderr, "                          train:    Do only training.  Input should\n");
  fprintf(stderr, "                                    be multiple FASTA of one or more\n");
  fprintf(stderr, "                                    closely related genomes.  Output\n");
  fprintf(stderr, "                                    is a training file.\n");
  fprintf(stderr, "                          anon:     Anonymous sequences, analyze using\n");
  fprintf(stderr, "                                    preset training files, ideal for\n");
  fprintf(stderr, "                                    metagenomic data or single short\n");
  fprintf(stderr, "                                    sequences.\n");
  fprintf(stderr, "  -g, --trans_table:    Specify a translation table to use\n");
  fprintf(stderr, "                          0:    Auto-detect (Default)\n");
  fprintf(stderr, "                          11:   Standard Bacterial/Archaeal\n");
  fprintf(stderr, "                          4:    Mycoplasma/Spiroplasma\n");
  fprintf(stderr, "                          #:    Other genetic codes 1-25\n");
  fprintf(stderr, "  -c, --nopartial:      Closed ends.  Do not allow partial genes\n");
  fprintf(stderr, "                        (genes that run off edges or into gaps.)\n");
  fprintf(stderr, "  -n, --nogaps:         Do not treat runs of N's as gaps.  This option\n");
  fprintf(stderr, "                        will build gene models that span runs of N's.\n");
  fprintf(stderr, "  -t, --training_file:  Read and use the specified training file\n");
  fprintf(stderr, "                        instead of training on the input sequence(s)\n");
  fprintf(stderr, "                        (Only usable in normal mode.)\n");
  fprintf(stderr, "\nInput/Output Parameters\n\n");
  fprintf(stderr, "  -i, --input_file:     Specify input file (default stdin).\n");
  fprintf(stderr, "  -o, --output_file:    Specify output file (default stdout).\n");
  fprintf(stderr, "  -a, --protein_file:   Write protein translations to the named file.\n");
  fprintf(stderr, "  -d, --mrna_file:      Write nucleotide sequences of genes to the\n");
  fprintf(stderr, "                        named file.\n");
  fprintf(stderr, "  -s, --start_file:     Write all potential genes (with scores) to the\n");
  fprintf(stderr, "                        named file.\n");
  fprintf(stderr, "  -f, --output_format:  Specify output format (gbk, gff, sqn, or sco).\n");
  fprintf(stderr, "                          gff:  GFF format (Default)\n");
  fprintf(stderr, "                          gbk:  Genbank-like format\n");
  fprintf(stderr, "                          sqn:  Sequin feature table format\n");
  fprintf(stderr, "                          sco:  Simple coordinate output\n");
  fprintf(stderr, "  -q, --quiet:          Run quietly (suppress normal stderr output).\n");
  fprintf(stderr, "\nOther Parameters\n\n");
  fprintf(stderr, "  -h, --help:     Print help menu and exit.\n");
  fprintf(stderr, "  -v, --version:  Print version number and exit.\n\n");
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
