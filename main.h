/******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tenum_nodesessee / UT-Battelle

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

#ifndef _MAIN_H_
#define _MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "anonymous.h"
#include "datatypes.h"
#include "dprog.h"
#include "gene.h"
#include "node.h"
#include "training.h"

#define VERSION "3.0.0-rc.1"
#define DATE "February, 2016"
#define TEXTSIZE 10000

struct _option
{
  char letter;                 /* Single letter identifier for option */
  char optarg[TEXTSIZE];       /* Parameter (if option takes one) */
  int index;                   /* Index of next argument to parse */
};

void version();
void usage(char *);
void help();
int allocate_memory(unsigned char **, unsigned char **, unsigned char **,
                    struct _node **, struct _gene **, struct _gene_data **,
                    struct _preset_genome_bin *, struct _gene ***);
void parse_arguments(int, char **, char *, char *, char *, char *, char *,
                     char *, char *, int *, int *, int *, int *, int *,
                     int *, int *);
void get_option(int, char **, struct _option *);
void header(int, int);
void log_text(int, char *);
int detect_input_and_handle_windows_stdin(int, int, char *);
int copy_standard_input_to_file(char *, int);
void open_files(char *, char *, char *, char *, char *, char *, FILE **,
                FILE **, FILE **, FILE **, FILE **, FILE **);

void free_variables(unsigned char *, unsigned char *, unsigned char *,
                    struct _node *, struct _gene *,
                    struct _preset_genome_bin *, struct _gene **);
void close_filehandles(FILE *, FILE *, FILE *, FILE *, FILE *, FILE *);

#endif
