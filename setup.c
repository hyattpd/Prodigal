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

#include "setup.h"

/* Print version number and exit */
void version() {
  printf("\nProdigal v%s: %s\n\n", VERSION, DATE);
  exit(0);
}

/* Print usage information and exit */
void usage(char *msg) {
  fprintf(stderr, "\nError: %s\n", msg);
  fprintf(stderr, "\nUsage:  prodigal [-a protein_file] [-c] [-d mrna_file]");
  fprintf(stderr, " [-f out_format]\n");
  fprintf(stderr, "                 [-g trans_table] [-h] [-i input_file]");
  fprintf(stderr, " [-m mode] [-n]\n");
  fprintf(stderr, "                 [-o output_file] [-q] [-s start_file]");
  fprintf(stderr, " [-t train_file]\n                 [-v] [-w summ_file]\n");
  fprintf(stderr, "\nDo 'prodigal -h' for more information.\n\n");
  exit(2);
}

/* Print help message and exit */
void help() {
  printf("\nUsage:  prodigal [-a protein_file] [-c] [-d mrna_file]");
  printf(" [-f out_format]\n");
  printf("                 [-g trans_table] [-h] [-i input_file]");
  printf(" [-m mode] [-n]\n");
  printf("                 [-o output_file] [-q] [-s start_file]");
  printf(" [-t train_file]\n                 [-v] [-w summ_file]\n");
  printf("\nGene Modeling Parameters\n\n");
  printf("  -m, --mode:           Specify mode (normal, train, or anon).\n");
  printf("                          normal:   Single genome, any number of\n");
  printf("                                    sequences. (Default)\n");
  printf("                          train:    Do only training.  Input should\n");
  printf("                                    be multiple FASTA of one or more\n");
  printf("                                    closely related genomes.  Output\n");
  printf("                                    is a training file.\n");
  printf("                          anon:     Anonymous sequences, analyze using\n");
  printf("                                    preset training files, ideal for\n");
  printf("                                    metagenomic data or single short\n");
  printf("                                    sequences.\n");
  printf("  -g, --trans_table:    Specify a translation table to use\n");
  printf("                          Auto: Tries 11 then 4 (Default)\n");
  printf("                          11:   Standard Bacteria/Archaea\n");
  printf("                          4:    Mycoplasma/Spiroplasma\n");
  printf("                          #:    Other genetic codes 1-25\n");
  printf("  -c, --nopartial:      Closed ends.  Do not allow partial genes\n");
  printf("                        (genes that run off edges or into gaps.)\n");
  printf("  -n, --nogaps:         Do not treat runs of N's as gaps.  This option\n");
  printf("                        will build gene models that span runs of N's.\n");
  printf("  -t, --training_file:  Read and use the specified training file\n");
  printf("                        instead of training on the input sequence(s)\n");
  printf("                        (Only usable in normal mode.)\n");
  printf("\nInput/Output Parameters\n\n");
  printf("  -i, --input_file:     Specify input file (default stdin).\n");
  printf("  -o, --output_file:    Specify output file (default stdout).\n");
  printf("  -a, --protein_file:   Write protein translations to the named file.\n");
  printf("  -d, --mrna_file:      Write nucleotide sequences of genes to the\n");
  printf("                        named file.\n");
  printf("  -w, --summ_file:      Write summary statistics to the named file.\n");
  printf("  -s, --start_file:     Write all potential genes (with scores) to the\n");
  printf("                        named file.\n");
  printf("  -f, --output_format:  Specify output format (gbk, gff, sqn, or sco).\n");
  printf("                          gff:  GFF format (Default)\n");
  printf("                          gbk:  Genbank-like format\n");
  printf("                          sqn:  Sequin feature table format\n");
  printf("                          sco:  Simple coordinate output\n");
  printf("  -q, --quiet:          Run quietly (suppress normal stderr output).\n");
  printf("\nOther Parameters\n\n");
  printf("  -h, --help:     Print help menu and exit.\n");
  printf("  -v, --version:  Print version number and exit.\n\n");
  exit(0);
}

/* Allocates memory for data structures and memsets them all to 0 */
int initialize_data_structures(unsigned char **seq, unsigned char **rseq,
                                unsigned char **useq, struct _node **nodes,
                                struct _gene **genes, struct _training *tinf,
                                struct _preset_genome_bin *presets, struct
                                _summary *statistics) {
  int i;

  *seq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char));
  *rseq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char));
  *useq = (unsigned char *)malloc(MAX_SEQ/8*sizeof(unsigned char));
  *nodes = (struct _node *)malloc(STT_NOD*sizeof(struct _node));
  *genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene));
  if(*seq == NULL || *rseq == NULL || *nodes == NULL || *genes == NULL) 
    return -1;
  memset(*seq, 0, MAX_SEQ/4*sizeof(unsigned char));
  memset(*rseq, 0, MAX_SEQ/4*sizeof(unsigned char));
  memset(*useq, 0, MAX_SEQ/8*sizeof(unsigned char));
  memset(*nodes, 0, STT_NOD*sizeof(struct _node));
  memset(*genes, 0, MAX_GENES*sizeof(struct _gene));
  memset(tinf, 0, sizeof(struct _training));
  memset(statistics, 0, sizeof(struct _summary));

  for(i = 0; i < NUM_PRESET_GENOME; i++) {
    memset(&presets[i], 0, sizeof(struct _preset_genome_bin));
    strcpy(presets[i].desc, "None");
    presets[i].tinf = (struct _training *)malloc(sizeof(struct _training));
    if(presets[i].tinf == NULL) return -1;
    memset(presets[i].tinf, 0, sizeof(struct _training));
  }
  return 0;
}

/* Initialize argument variables, parse command line arguments, */
/* and validate the arguments for consistency. */
void parse_arguments(int argc, char **argv, char *input_file, char 
                     *output_file, char *train_file, char *amino_file, char
                     *nuc_file, char *start_file, char *summ_file, int *mode,
                     int *outfmt, int *genetic_code, int *closed, int 
                     *cross_gaps, int *quiet) {
  int i, j;

  input_file[0] = '\0'; output_file[0] = '\0'; train_file[0] = '\0';
  start_file[0] = '\0'; summ_file[0] = '\0'; 
  amino_file[0] = '\0'; nuc_file[0] = '\0'; 
  *mode = 0; *outfmt = -1; *genetic_code = -1; 
  *closed = 0; *cross_gaps = 0; *quiet = 0;

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
       strcmp(argv[i], "-w") == 0 || 
       strcmp(argv[i], "--summ_file") == 0 ||
       strcmp(argv[i], "-i") == 0 || 
       strcmp(argv[i], "--input_file") == 0 ||
       strcmp(argv[i], "-o") == 0 || 
       strcmp(argv[i], "--output_file") == 0 ||
       strcmp(argv[i], "-m") == 0 || 
       strcmp(argv[i], "-p") == 0 || 
       strcmp(argv[i], "--mode") == 0))
      usage("-a/-d/-f/-g/-i/-m/-o/-s/-t/-w options require valid parameters.");
    else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--nopartial") == 0)
      *closed = 1;
    else if(strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--quiet") == 0)
      *quiet = 1;
    else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--nogaps") == 0)
      *cross_gaps = 1;
    else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) 
      help();
    else if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) 
      version();
    else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--protein_file") 
            == 0) {
      strcpy(amino_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--mrna_file") == 0) {
      strcpy(nuc_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-i") == 0 || 
            strcmp(argv[i], "--input_file") == 0) {
      strcpy(input_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || 
            strcmp(argv[i], "--output_file") == 0) {
      strcpy(output_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || 
            strcmp(argv[i], "--start_file") == 0) {
      strcpy(start_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-w") == 0 || 
            strcmp(argv[i], "--summ_file") == 0) {
      strcpy(summ_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || 
            strcmp(argv[i], "--training_file") == 0) {
      strcpy(train_file, argv[i+1]);
      i++;
    }
    else if(strcmp(argv[i], "-g") == 0 || 
            strcmp(argv[i], "--trans_table") == 0) {
      *genetic_code = atoi(argv[i+1]);
      if(*genetic_code < 0 || *genetic_code > 25 || *genetic_code == 7
         || *genetic_code == 8 || (*genetic_code >= 17 && *genetic_code
         <= 20))
        usage("Invalid or unsupported genetic code specified.");
      i++;
    }
    else if(strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mode") == 0) {
      if(argv[i+1][0] == 'n') *mode = 0;
      else if(argv[i+1][0] == 't') *mode = 1; 
      else if(argv[i+1][0] == 'a') *mode = 2; 
      else usage("Invalid mode specified (should be normal, train, or anon).");
      i++;
    }
    else if(strcmp(argv[i], "-p") == 0) { /* deprecated but preserved atm */
      if(argv[i+1][0] == 's') *mode = 0;
      else if(argv[i+1][0] == 'm') *mode = 2; 
      else usage("Invalid procedure specified (should be single or meta).");
      fprintf(stderr, "Warning: '-p meta' is deprecated.  Should use ");
      fprintf(stderr, "'-m anon' instead for metagenomic sequences.\n");
      i++;
    }
    else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--output_format") 
            == 0) {
      if(strncmp(argv[i+1], "0", 1) == 0 || strcmp(argv[i+1], "gbk") == 0)
        *outfmt = 0;
      else if(strncmp(argv[i+1], "1", 1) == 0 || strcmp(argv[i+1], "gca") == 0)
        *outfmt = 1;
      else if(strncmp(argv[i+1], "2", 1) == 0 || strcmp(argv[i+1], "sco") == 0)
        *outfmt = 2;
      else if(strncmp(argv[i+1], "3", 1) == 0 || strcmp(argv[i+1], "gff") == 0)
        *outfmt = 3;
      else if(strncmp(argv[i+1], "4", 1) == 0 || strcmp(argv[i+1], "sqn") == 0)
        *outfmt = 4;
      else usage("Invalid output format specified.");
      i++;
    }
    else {
      fprintf(stderr, "Unknown option '%s'.", argv[i]);
      usage("");
    }
  }

  /* Validation of arguments checking for conflicting options */

  /* Training mode can't have output format or extra files specified. */
  if(*mode == 1 && (strlen(start_file) > 0 || strlen(nuc_file) > 0 ||
     strlen(amino_file) > 0 || strlen(summ_file) > 0 || *outfmt != -1)) {
    usage("-a/-d/-f/-s/-w options cannot be used in training mode.");
  }

  /* Normal/anonymous can't have training files specified. */
  if(*mode != 0 && strlen(train_file) > 0) {
    usage("Can only specify training file in normal mode.");
  }

  /* Anonymous mode can't have a specified value for genetic code. */
  /* Nor can normal mode if using a training file. */
  if((*mode == 2 || (*mode == 0 && strlen(train_file) > 0)) &&
     *genetic_code != -1) {
    usage("Can't specify translation table with anonymous mode or a training file.");
  }
}

/* Print the header */
void header(int mode) {
  fprintf(stderr, "-------------------------------------\n");
  fprintf(stderr, "PRODIGAL v%s [%s]         \n", VERSION, DATE);
  fprintf(stderr, "Univ of Tenn / Oak Ridge National Lab\n");
  fprintf(stderr, "Doug Hyatt, Loren Hauser, et al.     \n");
  fprintf(stderr, "-------------------------------------\n");
  if(mode == 0) fprintf(stderr, "Mode: Normal, Phase: Training\n");
  else if(mode == 1) fprintf(stderr, "Mode: Training, Phase: Training\n");
  else if(mode == 2) fprintf(stderr, "Mode: Anonymous, Phase: Training\n");
}

/* If we're in normal mode and not reading from a training file, */
/* then we have to make two passes over the sequence.  Since we  */
/* rewind after the first pass (something Windows can't do to    */
/* stdin), we copy stdin to a temp file so that we can rewind it */
/* in Windows. If there's nothing present on stdin, we print     */
/* the help message and exit. Returns a 1 if piped input is      */
/* detected, 0 otherwise.                                        */
int detect_input_and_handle_windows_stdin(int argc, int quiet, 
                                           char *input_file) {
  int fnum, piped = 0;
  char input_copy[MAX_LINE];
  struct stat fbuf;
  pid_t pid;

  pid = getpid();
  sprintf(input_copy, "tmp.prodigal.stdin.%d", pid);

  fnum = fileno(stdin);
  if(fstat(fnum, &fbuf) == -1) {
    fprintf(stderr, "\nError: can't fstat standard input.\n\n");
    exit(3);
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
    strcpy(input_file, input_copy);
  }
  return piped;
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

/* Open files and set file pointers.  Exit if any files throw an error. */
void open_files(char *input_file, char *output_file, char *start_file,
                char *amino_file, char *nuc_file, char *summ_file,
                FILE **input_ptr, FILE **output_ptr, FILE **start_ptr,
                FILE **amino_ptr, FILE **nuc_ptr, FILE **summ_ptr) {
  if(input_file[0] != '\0') {
    *input_ptr = fopen(input_file, "r");
    if(*input_ptr == NULL) { 
      perror("\nError: can't open input file."); exit(7);
    }
  }
  if(output_file[0] != '\0') {
    *output_ptr = fopen(output_file, "w");
    if(*output_ptr == NULL) {
      perror("\nError: can't open output file."); exit(8);
    }
  }
  if(start_file[0] != '\0') {
    *start_ptr = fopen(start_file, "w");
    if(*start_ptr == NULL) {
      perror("\nError: can't open start file."); exit(8);
    }
  }
  if(amino_file[0] != '\0') {
    *amino_ptr = fopen(amino_file, "w");
    if(*amino_ptr == NULL) {
      perror("\nError: can't open translation file."); exit(8);
    }
  }
  if(nuc_file[0] != '\0') {
    *nuc_ptr = fopen(nuc_file, "w");
    if(*nuc_ptr == NULL) {
      perror("\nError: can't open gene nucleotide file."); exit(8);
    }
  }
  if(summ_file[0] != '\0') {
    *summ_ptr = fopen(summ_file, "w");
    if(*summ_ptr == NULL) {
      perror("\nError: can't open summary file."); exit(8);
    }
  }
}
