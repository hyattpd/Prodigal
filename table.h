/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    table.h: routines and data structures for parsing NCBI-style 64-character
    genetic code strings

    Code Author:  Mingye Wang

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

#ifndef _TABLE_H
#define _TABLE_H

/* 1.1 P_*: trinuc()/is_a() order: {a,g,c,t}aa, {a,g,c,t}ag, ... */
/* Extracted here because useful elsewhere! */
#define P_T 3
#define P_G 1
#define P_C 2
#define P_A 0
#define P_3(x, y, z) ((z) << 4 | (y) << 2 | (x))
#define P_TAG P_3(P_T, P_A, P_G)
#define P_TAA P_3(P_T, P_A, P_A)
#define P_TGA P_3(P_T, P_G, P_A)
/* Conversion from and to {ncbieaa,sncbieaa} */
void eaa_to_table(char[64], const char[65], const char[65]);
void table_to_eaa(const char[64], char[65], char[65]);

/* Conversion from and to a predefined table number */
int id_to_table(char[64], int);
int table_to_id(const char[64]);

/* Variant for handling user input, which is potentially malformed */
int cmdline_eaa_to_table(char[64], const char[130]);

#include "table-data.hh"
extern const char predefined_tables[MAXTABLE+1][2][65];
#endif
