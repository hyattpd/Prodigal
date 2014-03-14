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
#include "node.h"
#include "gene.h"
#include "anonymous.h"
#include "sequence.h"
#include "dprog.h"

#define VERSION "2.7.0"
#define DATE "March, 2014"

#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000

void version();
void usage(char *);
void help();
int copy_standard_input_to_file(char *, int);
int initialize_data_structures(unsigned char **, unsigned char **, unsigned 
                               char **, struct _node **, struct _gene **, 
                               struct _training *, struct 
                               _preset_genome_bin *);
void parse_arguments(int, char **, char **, char **, char **, char **, char **,
                     char **, char **, int *, int *, int *, int *, int *,
                     int *);

#endif
