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

#ifndef SETUP_H_
#define SETUP_H_

#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "datatypes.h"
#include "summary.h"
#include "gene.h"
#include "anonymous.h"

#define VERSION "2.7.0rc1"
#define DATE "May, 2014"

void version();
void usage(char *);
void help();

int initialize_data_structures(unsigned char **, unsigned char **, unsigned 
                               char **, struct _node **, struct _gene **, 
                               struct _training *, struct 
                               _preset_genome_bin *, struct _summary *);
void parse_arguments(int, char **, char *, char *, char *, char *, char *,
                     char *, char *, int *, int *, int *, int *, int *,
                     int *);
void header(int);
int detect_input_and_handle_windows_stdin(int, int, char *);
int copy_standard_input_to_file(char *, int);
void open_files(char *, char *, char *, char *, char *, char *, FILE **,
                FILE **, FILE **, FILE **, FILE **, FILE **);

#endif
