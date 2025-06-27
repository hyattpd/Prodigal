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

#include "table.h"
#include <string.h>

/* Part 1. Convert the NCBI amino-acid order from/to trinuc() order. */
/* 1.2 N_*: ncbi codon order: tt{t,c,a,g}, tc{t,c,a,g}, ... */
#define N_T 0
#define N_C 1
#define N_A 2
#define N_G 3
#define N_3(x, y, z) ((x) << 4 | (y) << 2 | (z))

const char ncbi_to_trinuc[64] = {
  P_3(P_T, P_T, P_T), P_3(P_T, P_T, P_C), P_3(P_T, P_T, P_A), P_3(P_T, P_T, P_G),
  P_3(P_T, P_C, P_T), P_3(P_T, P_C, P_C), P_3(P_T, P_C, P_A), P_3(P_T, P_C, P_G),
  P_3(P_T, P_A, P_T), P_3(P_T, P_A, P_C), P_3(P_T, P_A, P_A), P_3(P_T, P_A, P_G),
  P_3(P_T, P_G, P_T), P_3(P_T, P_G, P_C), P_3(P_T, P_G, P_A), P_3(P_T, P_G, P_G),
  P_3(P_C, P_T, P_T), P_3(P_C, P_T, P_C), P_3(P_C, P_T, P_A), P_3(P_C, P_T, P_G),
  P_3(P_C, P_C, P_T), P_3(P_C, P_C, P_C), P_3(P_C, P_C, P_A), P_3(P_C, P_C, P_G),
  P_3(P_C, P_A, P_T), P_3(P_C, P_A, P_C), P_3(P_C, P_A, P_A), P_3(P_C, P_A, P_G),
  P_3(P_C, P_G, P_T), P_3(P_C, P_G, P_C), P_3(P_C, P_G, P_A), P_3(P_C, P_G, P_G),
  P_3(P_A, P_T, P_T), P_3(P_A, P_T, P_C), P_3(P_A, P_T, P_A), P_3(P_A, P_T, P_G),
  P_3(P_A, P_C, P_T), P_3(P_A, P_C, P_C), P_3(P_A, P_C, P_A), P_3(P_A, P_C, P_G),
  P_3(P_A, P_A, P_T), P_3(P_A, P_A, P_C), P_3(P_A, P_A, P_A), P_3(P_A, P_A, P_G),
  P_3(P_A, P_G, P_T), P_3(P_A, P_G, P_C), P_3(P_A, P_G, P_A), P_3(P_A, P_G, P_G),
  P_3(P_G, P_T, P_T), P_3(P_G, P_T, P_C), P_3(P_G, P_T, P_A), P_3(P_G, P_T, P_G),
  P_3(P_G, P_C, P_T), P_3(P_G, P_C, P_C), P_3(P_G, P_C, P_A), P_3(P_G, P_C, P_G),
  P_3(P_G, P_A, P_T), P_3(P_G, P_A, P_C), P_3(P_G, P_A, P_A), P_3(P_G, P_A, P_G),
  P_3(P_G, P_G, P_T), P_3(P_G, P_G, P_C), P_3(P_G, P_G, P_A), P_3(P_G, P_G, P_G)
};

/* 1.1 P_*: trinuc()/is_a() order: {a,g,c,t}aa, {a,g,c,t}ga, ... */
const char trinuc_to_ncbi[64] = {
  N_3(N_A, N_A, N_A), N_3(N_G, N_A, N_A), N_3(N_C, N_A, N_A), N_3(N_T, N_A, N_A),
  N_3(N_A, N_G, N_A), N_3(N_G, N_G, N_A), N_3(N_C, N_G, N_A), N_3(N_T, N_G, N_A),
  N_3(N_A, N_C, N_A), N_3(N_G, N_C, N_A), N_3(N_C, N_C, N_A), N_3(N_T, N_C, N_A),
  N_3(N_A, N_T, N_A), N_3(N_G, N_T, N_A), N_3(N_C, N_T, N_A), N_3(N_T, N_T, N_A),
  N_3(N_A, N_A, N_G), N_3(N_G, N_A, N_G), N_3(N_C, N_A, N_G), N_3(N_T, N_A, N_G),
  N_3(N_A, N_G, N_G), N_3(N_G, N_G, N_G), N_3(N_C, N_G, N_G), N_3(N_T, N_G, N_G),
  N_3(N_A, N_C, N_G), N_3(N_G, N_C, N_G), N_3(N_C, N_C, N_G), N_3(N_T, N_C, N_G),
  N_3(N_A, N_T, N_G), N_3(N_G, N_T, N_G), N_3(N_C, N_T, N_G), N_3(N_T, N_T, N_G),
  N_3(N_A, N_A, N_C), N_3(N_G, N_A, N_C), N_3(N_C, N_A, N_C), N_3(N_T, N_A, N_C),
  N_3(N_A, N_G, N_C), N_3(N_G, N_G, N_C), N_3(N_C, N_G, N_C), N_3(N_T, N_G, N_C),
  N_3(N_A, N_C, N_C), N_3(N_G, N_C, N_C), N_3(N_C, N_C, N_C), N_3(N_T, N_C, N_C),
  N_3(N_A, N_T, N_C), N_3(N_G, N_T, N_C), N_3(N_C, N_T, N_C), N_3(N_T, N_T, N_C),
  N_3(N_A, N_A, N_T), N_3(N_G, N_A, N_T), N_3(N_C, N_A, N_T), N_3(N_T, N_A, N_T),
  N_3(N_A, N_G, N_T), N_3(N_G, N_G, N_T), N_3(N_C, N_G, N_T), N_3(N_T, N_G, N_T),
  N_3(N_A, N_C, N_T), N_3(N_G, N_C, N_T), N_3(N_C, N_C, N_T), N_3(N_T, N_C, N_T),
  N_3(N_A, N_T, N_T), N_3(N_G, N_T, N_T), N_3(N_C, N_T, N_T), N_3(N_T, N_T, N_T)
};

/* We only do ATG, GTG, TTG start, and it's only ever Met */
const char *sncbieaa_template = 
  "---M---------------M---------------M----------------------------";
inline static char sncbieaa_get(const char sncbieaa[65], int i) {
  if (sncbieaa[i] != '*' && sncbieaa[i] != '-')
    return sncbieaa_template[i];
  else
    return sncbieaa[i];
}

void eaa_to_table(char out[64], const char ncbieaa[65], const char sncbieaa[65]) {
  char temp_table[64];
  memcpy(temp_table, ncbieaa, 64);
  for (int i = 0; i < 64; i++) {
    /* The content of the table is the same as ncbieaa, except we encode the
       data from sncbieaa in bits. */
    if (sncbieaa[i] == '*')    /* Make all conditional stops into simple '*' */
        temp_table[i] = '*';   /* We could encode it as 0x80 but then what? */
    else if (sncbieaa[i] != '-')
        temp_table[i] |= 0x20; 

    out[(unsigned) ncbi_to_trinuc[i]] = temp_table[i];
  }
}

void table_to_eaa(const char table[64], char ncbieaa[65], char sncbieaa[65]) {
  for (int i = 0; i < 64; i++) {
    unsigned trinuc = trinuc_to_ncbi[i];
    char amino = table[trinuc];
    ncbieaa[i] = amino &~ 0xA0;
    if (amino & 0x80 || ncbieaa[i] == '*')
      sncbieaa[i] = '*';
    else if (amino & 0x20)
      sncbieaa[i] = 'M';
    else
      sncbieaa[i] = '-';
  }
  ncbieaa[64] = '\0';
  sncbieaa[64] = '\0';
}

/* Part 2. Predefined tables, which is a mix of NCBI and other data.
    */
#include "table-data.cc"

int id_to_table(char out[64], int id) {
  if (id < 0 || id > MAXTABLE) {
    return -1;
  }
  eaa_to_table(out, predefined_tables[id][0], predefined_tables[id][1]);
  return 0;
}

int table_to_id(const char table[64]) {
  char ncbieaa[65], sncbieaa[65];
  table_to_eaa(table, ncbieaa, sncbieaa);
  for (int i = 0; i <= MAXTABLE; i++) {
    if (strcmp(ncbieaa, predefined_tables[i][0]) == 0 &&
        strcmp(ncbieaa, predefined_tables[i][1]) == 0) {
      return i;
    }
  }
  return -1;
}

/* Part 3. Handling user input */
/* Unfortunately * also has 0x20 set, so we use the unset version... */
const char *alphabet = "ACDEFGHIKLMNPQRSTVWYUO\x0A";
const char *salphabet = "ACDEFGHIKLMNPQRSTVWYUO*-";
int cmdline_eaa_to_table(char table[64], const char input[130]) {
  int len = strlen(input);
  int input_is_table = 0;
  char ncbieaa[65], sncbieaa[65];
  const char *inputs = NULL;
  if (len != 64 && len != 129) {
    return -1;
  }

  for (int i = 0; i < 64; i++) {
    if (strchr(alphabet, (input[i] & (~0xA0))) == NULL)
      return -2;
    if (input[i] != '*' && (input[i] & 0xA0) != 0)
      input_is_table = 1;
  }
  if (input_is_table) {
    memcpy(table, input, 64);
    return 0;
  }

  memcpy(ncbieaa, input, 64);
  ncbieaa[64] = '\0';
  if (len == 64) {
    memcpy(sncbieaa, sncbieaa_template, 64);
    sncbieaa[64] = '\0';
    for (int i = 0; i < 64; i++) {
      if (ncbieaa[i] == '*') {
        sncbieaa[i] = '*';
      }
    }
  } else {
    inputs = input + 64;
    for (int i = 0; i < 64; i++) {
      sncbieaa[i] = sncbieaa_get(inputs, i);
      if (strchr(salphabet, inputs[i]) == NULL)
        return -3;
    }
    sncbieaa[64] = '\0';
  }

  eaa_to_table(table, ncbieaa, sncbieaa);
  return 0;
}

#ifdef TABLE_UTIL
#include <stdio.h>
#include <stdlib.h>

void usage(char* a0) {
  fprintf(stderr, "Usage: %s ce2t NCBIEAA[,SNCBIEAA]  "
    "# parse -g NCBI\n", a0);
  fprintf(stderr, "       %s e2t  NCBIEAA SNCBIEAA    "
    "# parse ordinary NCBI\n", a0);
  fprintf(stderr, "       %s t2e  TABLE  # table string to NCBI\n", a0);
  fprintf(stderr, "       %s id2t ID     # id to table\n", a0);
  fprintf(stderr, "       %s t2id TABLE  # table to id\n\n", a0);
}
int main(int argc, char *argv[]) {
  if (argc < 3) {
    usage(argv[0]);
    return 1;
  }
  char table[65], ncbieaa[65], sncbieaa[65];
  int ret;
  if (strcmp(argv[1], "ce2t") == 0) {
    if (ret = cmdline_eaa_to_table(table, argv[2])) {
      fprintf(stderr, "cmdline_eaa_to_table: %d\n", ret);
      usage(argv[0]);
      return 1;
    }
    table_to_eaa(table, ncbieaa, sncbieaa);
    printf("%s,%s\n", ncbieaa, sncbieaa);
  } else if (strcmp(argv[1], "e2t") == 0) {
    eaa_to_table(table, argv[2], argv[3]);
    table[64] = '\0';
    puts(table);
  } else if (strcmp(argv[1], "t2e") == 0) {
    table_to_eaa(argv[2], ncbieaa, sncbieaa);
    printf("%s,%s\n", ncbieaa, sncbieaa);
  } else if (strcmp(argv[1], "id2t") == 0) {
    int id = atoi(argv[2]);
    if (ret = id_to_table(table, id)) {
      fprintf(stderr, "id_to_table: %d\n", ret);
      usage(argv[0]);
      return 1;
    }
    printf("%.*s\n", 64, table);
  } else if (strcmp(argv[1], "t2id") == 0) {
    int id = table_to_id(argv[2]);
    printf("%d\n", id);
  } else {
    usage(argv[0]);
    return 1;
  }
}
#endif
