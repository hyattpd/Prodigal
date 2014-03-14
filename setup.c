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

void version() {
  fprintf(stderr, "\nProdigal V%s: %s\n\n", VERSION, DATE);
  exit(0);
}

void usage(char *msg) {
  fprintf(stderr, "\nError: %s\n", msg);
  fprintf(stderr, "\nUsage:  prodigal [-a protein_file] [-c] [-d mrna_file]");
  fprintf(stderr, " [-f out_format]\n");
  fprintf(stderr, "                 [-g trans_table] [-h] [-i input_file]");
  fprintf(stderr, " [-m mode] [-n]\n");
  fprintf(stderr, "                 [-o output_file] [-q] [-s start_file]");
  fprintf(stderr, " [-t train_file]\n                 [-v] [-w summ_file]\n");
  fprintf(stderr, "\nDo 'prodigal -h' for more information.\n\n");
  exit(15);
}

void help() {
  fprintf(stderr, "\nUsage:  prodigal [-a protein_file] [-c] [-d mrna_file]");
  fprintf(stderr, " [-f out_format]\n");
  fprintf(stderr, "                 [-g trans_table] [-h] [-i input_file]");
  fprintf(stderr, " [-m mode] [-n]\n");
  fprintf(stderr, "                 [-o output_file] [-q] [-s start_file]");
  fprintf(stderr, " [-t train_file]\n                 [-v] [-w summ_file]\n");
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
  fprintf(stderr, "                          Auto: Tries 11 then 4 (Default)\n");
  fprintf(stderr, "                          11:   Standard Bacteria/Archaea\n");
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
  fprintf(stderr, "  -w, --summ_file:      Write summary statistics to the named file.\n");
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

/* Allocates memory for data structures and memsets them all to 0 */
int initialize_data_structures(unsigned char **seq, unsigned char **rseq,
                                unsigned char **useq, struct _node **nodes,
                                struct _gene **genes, struct _training *tinf,
                                struct _preset_genome_bin *presets) {
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

  for(i = 0; i < NUM_PRESET_GENOME; i++) {
    memset(&presets[i], 0, sizeof(struct _preset_genome_bin));
    strcpy(presets[i].desc, "None");
    presets[i].tinf = (struct _training *)malloc(sizeof(struct _training));
    if(presets[i].tinf == NULL) return -1;
    memset(presets[i].tinf, 0, sizeof(struct _training));
  }
  return 0;
}

/* Parse command line arguments */
void parse_arguments(int argc, char **argv, char **input_file, char 
                     **output_file, char **train_file, char **amino_file, char
                     **nuc_file, char **start_file, char **summ_file, int *mode,
                     int *outfmt, int *genetic_code, int *closed, int 
                     *cross_gaps, int *quiet) {
  int i, j;

  *input_file = NULL;
  *output_file = NULL;
  *train_file = NULL;
  *start_file = NULL;
  *amino_file = NULL;
  *nuc_file = NULL;
  *summ_file = NULL;
  *mode = 0;
  *outfmt = -1;
  *genetic_code = -1;
  *closed = 0;
  *cross_gaps = 0;
  *quiet = 0;

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
      *amino_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--mrna_file") == 0) {
      *nuc_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-i") == 0 || 
            strcmp(argv[i], "--input_file") == 0) {
      *input_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-o") == 0 || 
            strcmp(argv[i], "--output_file") == 0) {
      *output_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-s") == 0 || 
            strcmp(argv[i], "--start_file") == 0) {
      *start_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-w") == 0 || 
            strcmp(argv[i], "--summ_file") == 0) {
      *summ_file = argv[i+1];
      i++;
    }
    else if(strcmp(argv[i], "-t") == 0 || 
            strcmp(argv[i], "--training_file") == 0) {
      *train_file = argv[i+1];
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
      fprintf(stderr, "Warning: '-p meta' is deprecated.  Should use -m anon ");
      fprintf(stderr, "instead.\n");
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
  if(*mode == 1 && (*start_file != NULL || *nuc_file != NULL ||
     *amino_file != NULL || *summ_file != NULL || *outfmt != -1)) {
    usage("-a/-d/-f/-s/-w options cannot be used in training mode.");
  }

  /* Normal/anonymous can't have training files specified. */
  if(*mode != 0 && *train_file != NULL) {
    usage("Can only specify training file in normal mode.");
  }

  /* Anonymous mode can't have a specified value for genetic code. */
  /* Nor can normal mode if using a training file. */
  if((*mode == 2 || (*mode == 0 && *train_file != NULL)) &&
     *genetic_code != -1) {
    usage("Can't specify translation table with anonymous mode or a training file.");
  }
}
