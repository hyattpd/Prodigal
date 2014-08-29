/******************************************************************************
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
******************************************************************************/

#include "sequence.h"

/******************************************************************************
  Read the sequence for training purposes.  If we encounter multiple
  sequences, we insert gaps or stop codons in between each one to allow/forbid
  training on partial genes.  When we hit MAX_SEQ bp, we stop and return what
  we've got so far for training.  This routine reads in FASTA, and has a very
  'loose' Genbank and Embl parser, but, to be safe, FASTA should generally be
  preferred.
******************************************************************************/
int read_seq_training(FILE *fp, unsigned char *seq, unsigned char *useq,
                      double *gc, int *num_seq)
{
  int i = 0;
  int state = 0;                /* 0 if haven't seen header, 1 if have */
  int num_header = 0;           /* Number of headers we've seen */
  int bit_len = 0;              /* Length of bitstring (2x actual slen) */
  int len = 0;                  /* Length of sequence created */
  int warn = 0;                 /* If we've warned already, 1 or 0 */
  int gc_cont = 0;              /* GC content of sequence */
  int gap_size = 0;             /* Size of gap parsed from Gbk file */
  char line[MAX_LINE+1] = "";   /* Buffer to read the sequence into */

  line[MAX_LINE] = '\0';
  while (fgets(line, MAX_LINE, fp) != NULL)
  {
    if (state == 0 && line[strlen(line)-1] != '\n' && warn == 0)
    {
      warn = 1;
      fprintf(stderr, "\n\nWarning: saw non-sequence line longer than ");
      fprintf(stderr, "%d chars, sequence might not be read ", MAX_LINE);
      fprintf(stderr, "correctly.\n\n");
    }
    if (line[0] == '>' || (line[0] == 'S' && line[1] == 'Q') ||
        (strlen(line) > 6 && strncmp(line, "ORIGIN", 6) == 0))
    {
      state = 1;
      /* Put a gap in between each training sequence */
      if (num_header > 0)
      {
        for (i = 0; i < 8; i++)
        {
          set(useq, len);
          bit_len+=2;
          len++;
        }
      }
      num_header++;
    }
    else if (state == 1 && (line[0] == '/' && line[1] == '/'))
    {
      state = 0;
    }
    else if (state == 1)
    {
      /* Some genbank files have Expand gap directives to process */
      if (strstr(line, "Expand") != NULL && strstr(line, "gap") != NULL)
      {
        sscanf(strstr(line, "gap")+4, "%d", &gap_size);
        if (gap_size < 1 || gap_size > MAX_LINE)
        {
          fprintf(stderr, "Error: gap size in gbk file can't exceed line");
          fprintf(stderr, " size.\n");
          exit(14);
        }
        for (i = 0; i < gap_size; i++)
        {
          line[i] = 'n';
        }
        line[i] = '\0';
      }
      for (i = 0; i < (int)strlen(line); i++)
      {
        if (line[i] < 'A' || line[i] > 'z')
        {
          continue;
        }
        if (line[i] == 'g' || line[i] == 'G')
        {
          set(seq, bit_len);
          gc_cont++;
        }
        else if (line[i] == 't' || line[i] == 'T' ||
                 line[i] == 'u' || line[i] == 'U')
        {
          set(seq, bit_len);
          set(seq, bit_len+1);
        }
        else if (line[i] == 'c' || line[i] == 'C')
        {
          set(seq, bit_len+1);
          gc_cont++;
        }
        else if (line[i] != 'a' && line[i] != 'A')
        {
          set(seq, bit_len+1);
          set(useq, len);
        }
        bit_len+=2;
        len++;
      }
    }
    if (len+MAX_LINE >= MAX_SEQ)
    {
      fprintf(stderr, "\n\nWarning: Sequence is long (max %d for training).\n",
              MAX_SEQ);
      fprintf(stderr, "Training on the first %d bases.\n\n", MAX_SEQ);
      break;
    }
  }
  *gc = ((double)gc_cont / (double)len);

  /* Exit if there are errors, warn if sequence is small */
  if (len == 0)
  {
    fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
    fprintf(stderr, "Genbank, or EMBL format).\n\n");
    exit(9);
  }
  if (len < MIN_SINGLE_GENOME)
  {
    fprintf(stderr, "\n\nError: Sequence must be %d", MIN_SINGLE_GENOME);
    fprintf(stderr, " characters (only %d read).\n(Consider", len);
    fprintf(stderr, " running with the '-p anon' option or finding");
    fprintf(stderr, " more contigs from the same genome.)\n\n");
    exit(10);
  }
  if (len < IDEAL_SINGLE_GENOME)
  {
    fprintf(stderr, "\n\nWarning: Ideally Prodigal should be given at");
    fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
    fprintf(stderr, "training.\nYou may get better results with the ");
    fprintf(stderr, "'-p anon' option.\n\n");
  }

  *num_seq = num_header;
  return len;
}

/* This routine reads in the next sequence in a FASTA/GB/EMBL file */
int next_seq_multi(FILE *fp, unsigned char *seq, unsigned char *useq,
                   int *seq_ctr, double *gc, char *cur_hdr, char *new_hdr)
{
  int i = 0;
  int reading_seq = 0;          /* State variable for if reading seq */
  int genbank_end = 0;          /* State variable for if at end of gbk */
  int bit_len = 0;              /* Length of bitstring (2x actual slen) */
  int len = 0;                  /* Length of sequence created */
  int warn = 0;                 /* If we've warned already, 1 or 0 */
  int gc_cont = 0;              /* GC content of sequence */
  int gap_size = 0;             /* Size of gap parsed from Gbk file */
  char line[MAX_LINE+1] = "";   /* Buffer to read sequence into */

  sprintf(new_hdr, "Prodigal_Seq_%d", *seq_ctr+2);
  if (*seq_ctr > 0)
  {
    reading_seq = 1;
  }
  line[MAX_LINE] = '\0';
  while (fgets(line, MAX_LINE, fp) != NULL)
  {
    if (reading_seq == 0 && line[strlen(line)-1] != '\n' && warn == 0)
    {
      warn = 1;
      fprintf(stderr, "\n\nWarning: saw non-sequence line longer than ");
      fprintf(stderr, "%d chars, sequence might not be read ", MAX_LINE);
      fprintf(stderr, "correctly.\n\n");
    }
    if (strlen(line) > 10 && strncmp(line, "DEFINITION", 10) == 0)
    {
      if (genbank_end == 0)
      {
        strcpy(cur_hdr, line+12);
        cur_hdr[strlen(cur_hdr)-1] = '\0';
      }
      else
      {
        strcpy(new_hdr, line+12);
        new_hdr[strlen(new_hdr)-1] = '\0';
      }
    }
    if (line[0] == '>' || (line[0] == 'S' && line[1] == 'Q') ||
        (strlen(line) > 6 && strncmp(line, "ORIGIN", 6) == 0))
    {
      if (reading_seq == 1 || genbank_end == 1 || *seq_ctr > 0)
      {
        if (line[0] == '>')
        {
          strcpy(new_hdr, line+1);
          new_hdr[strlen(new_hdr)-1] = '\0';
        }
        break;
      }
      if (line[0] == '>')
      {
        strcpy(cur_hdr, line+1);
        cur_hdr[strlen(cur_hdr)-1] = '\0';
      }
      reading_seq = 1;
    }
    else if (reading_seq == 1 && (line[0] == '/' && line[1] == '/'))
    {
      reading_seq = 0;
      genbank_end = 1;
    }
    else if (reading_seq == 1)
    {
      /* Handle gap directives in genbank file */
      if (strstr(line, "Expand") != NULL && strstr(line, "gap") != NULL)
      {
        sscanf(strstr(line, "gap")+4, "%d", &gap_size);
        if (gap_size < 1 || gap_size > MAX_LINE)
        {
          fprintf(stderr, "Error: gap size in gbk file can't exceed line");
          fprintf(stderr, " size.\n");
          exit(15);
        }
        for (i = 0; i < gap_size; i++)
        {
          line[i] = 'n';
        }
        line[i] = '\0';
      }
      for (i = 0; i < (int)strlen(line); i++)
      {
        if (line[i] < 'A' || line[i] > 'z')
        {
          continue;
        }
        if (line[i] == 'g' || line[i] == 'G')
        {
          set(seq, bit_len);
          gc_cont++;
        }
        else if (line[i] == 't' || line[i] == 'T' ||
                 line[i] == 'u' || line[i] == 'U')
        {
          set(seq, bit_len);
          set(seq, bit_len+1);
        }
        else if (line[i] == 'c' || line[i] == 'C')
        {
          set(seq, bit_len+1);
          gc_cont++;
        }
        else if (line[i] != 'a' && line[i] != 'A')
        {
          set(seq, bit_len+1);
          set(useq, len);
        }
        bit_len+=2;
        len++;
      }
    }
    if (len+MAX_LINE >= MAX_SEQ)
    {
      fprintf(stderr, "Sequence too long (max %d permitted).\n", MAX_SEQ);
      exit(16);
    }
  }
  if (len == 0)
  {
    return -1;
  }
  *gc = ((double)gc_cont / (double)len);
  *seq_ctr = *seq_ctr + 1;
  return len;
}

/* Takes first word of header.  If the length of the header is 0, */
/* it instead constructs a header of the form Prodigal_Seq_<ID> */
void calc_short_header(char *header, char *short_header, int seq_ctr)
{
  int i = 0;

  strcpy(short_header, header);
  for (i = 0; i < (int)strlen(header); i++)
  {
    if (header[i] == ' ' || header[i] == '\t' || header[i] == '\r' ||
        header[i] == '\n')
    {
      strncpy(short_header, header, i);
      short_header[i] = '\0';
      break;
    }
  }
  if (i == 0)
  {
    sprintf(short_header, "Prodigal_Seq_%d", seq_ctr);
  }
}

/* Stores reverse complement of seq in rseq */
void reverse_seq(unsigned char *seq, unsigned char *rseq, unsigned char *useq,
                 int len)
{
  int i = 0;
  int bit_len = len*2;         /* Bitstring size is twice the seq length */
  for (i = 0; i < bit_len; i++)
  {
    if (test(seq, i) == 0)
    {
       set(rseq, bit_len-i-1+(i%2==0?-1:1));
    }
  }
  for (i = 0; i < len; i++)
  {
    if (test(useq, i) == 1)
    {
      toggle(rseq, bit_len-1-i*2);
      toggle(rseq, bit_len-2-i*2);
    }
  }
}

/* Simple routines to say whether or not bases are */
/* a, c, t, g, starts, stops, etc. */
int is_a(unsigned char *seq, int n)
{
  int index = n*2;
  if (test(seq, index) == 1 || test(seq, index+1) == 1)
  {
    return 0;
  }
  return 1;
}

int is_c(unsigned char *seq, int n)
{
  int index = n*2;
  if (test(seq, index) == 1 || test(seq, index+1) == 0)
  {
    return 0;
  }
  return 1;
}

int is_g(unsigned char *seq, int n)
{
  int index = n*2;
  if (test(seq, index) == 0 || test(seq, index+1) == 1)
  {
     return 0;
  }
  return 1;
}

int is_t(unsigned char *seq, int n)
{
  int index = n*2;
  if (test(seq, index) == 0 || test(seq, index+1) == 0)
  {
    return 0;
  }
  return 1;
}

int is_n(unsigned char *useq, int n)
{
  if (test(useq, n) == 0)
  {
    return 0;
  }
  return 1;
}

int is_stop(unsigned char *seq, int n, int trans_table)
{
  /* TAG: Not a stop in genetic codes 6, 15, 16, and 22 */
  if (is_t(seq, n) == 1 && is_a(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    if (trans_table == 6 || trans_table == 15 ||
        trans_table == 16 || trans_table == 22)
    {
       return 0;
    }
    return 1;
  }

  /* TGA: Not a stop in genetic codes 2-5, 9-10, 13-14, 21, 24-25 */
  if (is_t(seq, n) == 1 && is_g(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if ((trans_table >= 2 && trans_table <= 5) ||
        trans_table == 9 || trans_table == 10 ||
        trans_table == 13 || trans_table == 14 ||
        trans_table == 21 || trans_table == 24 ||
        trans_table == 25)
    {
      return 0;
    }
    return 1;
  }

  /* TAA: Not a stop in genetic codes 6 and 14 */
  if (is_t(seq, n) == 1 && is_a(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 6 || trans_table == 14)
    {
      return 0;
    }
    return 1;
  }

  /* Code 2: AGG and AGA are stops */
  if (trans_table == 2 && is_a(seq, n) == 1 && is_g(seq, n+1) == 1 &&
      is_a(seq, n+2) == 1)
  {
    return 1;
  }
  if (trans_table == 2 && is_a(seq, n) == 1 && is_g(seq, n+1) == 1 &&
      is_g(seq, n+2) == 1)
  {
    return 1;
  }

  /* Code 22: TCA is a stop */
  if (trans_table == 22 && is_t(seq, n) == 1 && is_c(seq, n+1) == 1 &&
      is_a(seq, n+2) == 1)
  {
    return 1;
  }

  /* Code 23: TTA is a stop */
  if (trans_table == 23 && is_t(seq, n) == 1 && is_t(seq, n+1) == 1 &&
      is_a(seq, n+2) == 1)
  {
    return 1;
  }

  return 0;
}

/* Prodigal only supports ATG/GTG/TTG starts as 'standard' possibilities. */
/* Some genetic codes have other initiation codons listed, but we do not  */
/* support these. */
int is_start(unsigned char *seq, int n, int trans_table)
{
  /* ATG: Always a start codon */
  if (is_a(seq, n) == 1 && is_t(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 1;
  }

  /* GTG: Start codon in 2/4/5/9/11/13/21/23/24/25 */
  if (is_g(seq, n) == 1 && is_t(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    if (trans_table == 2 || trans_table == 4 ||
        trans_table == 5 || trans_table == 9 ||
        trans_table == 11 || trans_table == 13 ||
        trans_table == 21 || trans_table == 23 ||
        trans_table == 24 || trans_table == 25)
    {
      return 1;
    }
  }

  /* TTG: Start codon in 4/5/11/13/24/25 */
  if (is_t(seq, n) == 1 && is_t(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    if (trans_table == 4 || trans_table == 5 ||
        trans_table == 11 || trans_table == 13 ||
        trans_table == 24 || trans_table == 25)
    {
      return 1;
    }
  }

  /* We do not handle other initiation codons */
  return 0;
}

/* Routines to say if codons are common start or stop triplets */
int is_atg(unsigned char *seq, int n)
{
  if (is_a(seq, n) == 0 || is_t(seq, n+1) == 0 || is_g(seq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

int is_gtg(unsigned char *seq, int n)
{
  if (is_g(seq, n) == 0 || is_t(seq, n+1) == 0 || is_g(seq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

int is_ttg(unsigned char *seq, int n)
{
  if (is_t(seq, n) == 0 || is_t(seq, n+1) == 0 || is_g(seq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

int is_taa(unsigned char *seq, int n)
{
  if (is_t(seq, n) == 0 || is_a(seq, n+1) == 0 || is_a(seq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

int is_tag(unsigned char *seq, int n)
{
  if (is_t(seq, n) == 0 || is_a(seq, n+1) == 0 || is_g(seq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

int is_tga(unsigned char *seq, int n)
{
  if (is_t(seq, n) == 0 || is_g(seq, n+1) == 0 || is_a(seq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

/* Returns 1 if entire codon is ambiguous */
int is_nnn(unsigned char *useq, int n)
{
  if (is_n(useq, n) == 0 || is_n(useq, n+1) == 0 || is_n(useq, n+2) == 0)
  {
    return 0;
  }
  return 1;
}

/* Returns 1 if any base in a codon is ambiguous */
int codon_has_n(unsigned char *useq, int n)
{
  if (is_n(useq, n) == 1 || is_n(useq, n+1) == 1 || is_n(useq, n+2) == 1)
  {
    return 1;
  }
  return 0;
}

/* Returns 1 if two entire codons to the left of this codon are ambiguous */
int gap_to_left(unsigned char *useq, int n)
{
  if (is_nnn(useq, n-3) == 1 && is_nnn(useq, n-6) == 1)
  {
    return 1;
  }
  else if (is_n(useq, n-3) == 1 && is_nnn(useq, n-6) == 1 &&
           is_nnn(useq, n-9) == 1)
  {
    return 1;
  }
  return 0;
}

/* Returns 1 if two entire codons to the right of this codon are ambiguous */
int gap_to_right(unsigned char *useq, int n)
{
  if (is_nnn(useq, n+3) == 1 && is_nnn(useq, n+6) == 1)
  {
    return 1;
  }
  else if (is_n(useq, n+5) == 1 && is_nnn(useq, n+6) == 1 &&
           is_nnn(useq, n+9) == 1)
  {
    return 1;
  }
  return 0;
}

/* Returns 1 if the base is a G or a C */
int is_gc(unsigned char *seq, int n)
{
  int index = n*2;
  if (test(seq, index) != test(seq, index+1))
  {
    return 1;
  }
  return 0;
}

/* Returns the probability a random codon should be a stop codon */
/* based on the GC content and genetic code of the organism */
double prob_stop(int tt, double gc)
{
  int i1 = 0;                     /* Loop variables, two for each base */
  int i2 = 0;                     /* (since manipulating bitstring */
  int i3 = 0;
  int i4 = 0;
  int i5 = 0;
  int i6 = 0;
  unsigned char codon[3] = "";    /* Current codon for consideration */
  double codon_prob = 0.0;        /* Prob of current codon */
  double stop_prob = 0.0;         /* Cumulative prob of any stop codon */

  for (i1 = 0; i1 < 6; i1++)
  {
    clear(codon, i1);
  }
  for (i1 = 0; i1 < 2; i1++)
  {
    if (i1 == 0)
    {
      clear(codon, 0);
    }
    else
    {
      set(codon, 0);
    }
    for (i2 = 0; i2 < 2; i2++)
    {
      if (i2 == 0)
      {
        clear(codon, 1);
      }
      else
      {
        set(codon, 1);
      }
      for (i3 = 0; i3 < 2; i3++)
      {
        if (i3 == 0)
        {
          clear(codon, 2);
        }
        else
        {
          set(codon, 2);
        }
        for (i4 = 0; i4 < 2; i4++)
        {
          if (i4 == 0)
          {
            clear(codon, 3);
          }
          else
          {
            set(codon, 3);
          }
          for (i5 = 0; i5 < 2; i5++)
          {
            if (i5 == 0)
            {
              clear(codon, 4);
            }
            else
            {
              set(codon, 4);
            }
            for (i6 = 0; i6 < 2; i6++)
            {
              if (i6 == 0)
              {
                clear(codon, 5);
              }
              else
              {
                set(codon, 5);
              }
              codon_prob = 1.0;
              if (is_gc(codon, 0) == 1)
              {
                codon_prob *= gc/2.0;
              }
              else
              {
                codon_prob *= (1.0-gc)/2.0;
              }
              if (is_gc(codon, 1) == 1)
              {
                codon_prob *= gc/2.0;
              }
              else
              {
                codon_prob *= (1.0-gc)/2.0;
              }
              if (is_gc(codon, 2) == 1)
              {
                codon_prob *= gc/2.0;
              }
              else
              {
                codon_prob *= (1.0-gc)/2.0;
              }
              if (is_stop(codon, 0, tt) == 1)
              {
                stop_prob += codon_prob;
              }
            }
          }
        }
      }
    }
  }
  return stop_prob;
}

/* Returns the GC content between 'n1' and 'n2' inclusive */
double gc_content(unsigned char *seq, int n1, int n2)
{
  int i = 0;
  double sum = 0.0;
  double gc = 0.0;
  for (i = n1; i <= n2; i++)
  {
    if (is_g(seq, i) == 1 || is_c(seq, i) == 1)
    {
      gc++;
    }
    sum++;
  }
  return gc/sum;
}

/* Returns a single amino acid for position 'n'. */
/* 'Is_init' being set indicates a start codon, which we always */
/* translate to 'M'. */
char amino(unsigned char *seq, int n, int trans_table, int is_init)
{
  if (is_stop(seq, n, trans_table) == 1)
  {
    return '*';
  }
  if (is_start(seq, n, trans_table) == 1 && is_init == 1)
  {
    return 'M';
  }
  if (is_t(seq, n) == 1 && is_t(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'F';
  }
  if (is_t(seq, n) == 1 && is_t(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'F';
  }
  if (is_t(seq, n) == 1 && is_t(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    return 'L';
  }
  if (is_t(seq, n) == 1 && is_t(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 'L';
  }
  if (is_t(seq, n) == 1 && is_c(seq, n+1) == 1)
  {
    return 'S';
  }
  if (is_t(seq, n) == 1 && is_a(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'Y';
  }
  if (is_t(seq, n) == 1 && is_a(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'Y';
  }
  if (is_t(seq, n) == 1 && is_a(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 6)
    {
      return 'Q';
    }
    if (trans_table == 14)
    {
      return 'Y';
    }
  }
  if (is_t(seq, n) == 1 && is_a(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    if (trans_table == 6 || trans_table == 15)
    {
      return 'Q';
    }
    if (trans_table == 16 || trans_table == 22)
    {
      return 'L';
    }
  }
  if (is_t(seq, n) == 1 && is_g(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'C';
  }
  if (is_t(seq, n) == 1 && is_g(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'C';
  }
  if (is_t(seq, n) == 1 && is_g(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 10)
    {
      return 'C';
    }
    if (trans_table == 25)
    {
      return 'G';
    }
    return 'W';
  }
  if (is_t(seq, n) == 1 && is_g(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 'W';
  }
  if (is_c(seq, n) == 1 && is_t(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    if (trans_table == 3)
    {
      return 'T';
    }
    return 'L';
  }
  if (is_c(seq, n) == 1 && is_t(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    if (trans_table == 3)
    {
      return 'T';
    }
    return 'L';
  }
  if (is_c(seq, n) == 1 && is_t(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 3)
    {
      return 'T';
    }
    return 'L';
  }
  if (is_c(seq, n) == 1 && is_t(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    if (trans_table == 3)
    {
      return 'T';
    }
    if (trans_table == 12)
    {
      return 'S';
    }
    return 'L';
  }
  if (is_c(seq, n) == 1 && is_c(seq, n+1) == 1)
  {
    return 'P';
  }
  if (is_c(seq, n) == 1 && is_a(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'H';
  }
  if (is_c(seq, n) == 1 && is_a(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'H';
  }
  if (is_c(seq, n) == 1 && is_a(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    return 'Q';
  }
  if (is_c(seq, n) == 1 && is_a(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 'Q';
  }
  if (is_c(seq, n) == 1 && is_g(seq, n+1) == 1)
  {
    return 'R';
  }
  if (is_a(seq, n) == 1 && is_t(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'I';
  }
  if (is_a(seq, n) == 1 && is_t(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'I';
  }
  if (is_a(seq, n) == 1 && is_t(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 2 || trans_table == 3 ||
        trans_table == 5 || trans_table == 13 ||
        trans_table == 21)
    {
      return 'M';
    }
    return 'I';
  }
  if (is_a(seq, n) == 1 && is_t(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 'M';
  }
  if (is_a(seq, n) == 1 && is_c(seq, n+1) == 1)
  {
    return 'T';
  }
  if (is_a(seq, n) == 1 && is_a(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'N';
  }
  if (is_a(seq, n) == 1 && is_a(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'N';
  }
  if (is_a(seq, n) == 1 && is_a(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 9 || trans_table == 14 ||
        trans_table == 21)
    {
      return 'N';
    }
    return 'K';
  }
  if (is_a(seq, n) == 1 && is_a(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 'K';
  }
  if (is_a(seq, n) == 1 && is_g(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'S';
  }
  if (is_a(seq, n) == 1 && is_g(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'S';
  }
  if (is_a(seq, n) == 1 && is_g(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    if (trans_table == 13)
    {
      return 'G';
    }
    if (trans_table == 5 || trans_table == 9 ||
        trans_table == 14 || trans_table == 21 ||
        trans_table == 24)
    {
      return 'S';
    }
    return 'R';
  }
  if (is_a(seq, n) == 1 && is_g(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    if (trans_table == 13)
    {
      return 'G';
    }
    if (trans_table == 5 || trans_table == 9 ||
        trans_table == 14 || trans_table == 21)
    {
      return 'S';
    }
    if (trans_table == 24)
    {
      return 'K';
    }
    return 'R';
  }
  if (is_g(seq, n) == 1 && is_t(seq, n+1) == 1)
  {
    return 'V';
  }
  if (is_g(seq, n) == 1 && is_c(seq, n+1) == 1)
  {
    return 'A';
  }
  if (is_g(seq, n) == 1 && is_a(seq, n+1) == 1 && is_t(seq, n+2) == 1)
  {
    return 'D';
  }
  if (is_g(seq, n) == 1 && is_a(seq, n+1) == 1 && is_c(seq, n+2) == 1)
  {
    return 'D';
  }
  if (is_g(seq, n) == 1 && is_a(seq, n+1) == 1 && is_a(seq, n+2) == 1)
  {
    return 'E';
  }
  if (is_g(seq, n) == 1 && is_a(seq, n+1) == 1 && is_g(seq, n+2) == 1)
  {
    return 'E';
  }
  if (is_g(seq, n) == 1 && is_g(seq, n+1) == 1)
  {
    return 'G';
  }
  return 'X';
}

/* Converts an amino acid letter to a numerical value */
int amino_num(char aa)
{
  if (aa == 'a' || aa == 'A')
  {
    return 0;
  }
  if (aa == 'c' || aa == 'C')
  {
    return 1;
  }
  if (aa == 'd' || aa == 'D')
  {
    return 2;
  }
  if (aa == 'e' || aa == 'E')
  {
    return 3;
  }
  if (aa == 'f' || aa == 'F')
  {
    return 4;
  }
  if (aa == 'g' || aa == 'G')
  {
    return 5;
  }
  if (aa == 'h' || aa == 'H')
  {
    return 6;
  }
  if (aa == 'i' || aa == 'I')
  {
    return 7;
  }
  if (aa == 'k' || aa == 'K')
  {
    return 8;
  }
  if (aa == 'l' || aa == 'L')
  {
    return 9;
  }
  if (aa == 'm' || aa == 'M')
  {
    return 10;
  }
  if (aa == 'n' || aa == 'N')
  {
    return 11;
  }
  if (aa == 'p' || aa == 'P')
  {
    return 12;
  }
  if (aa == 'q' || aa == 'Q')
  {
    return 13;
  }
  if (aa == 'r' || aa == 'R')
  {
    return 14;
  }
  if (aa == 's' || aa == 'S')
  {
    return 15;
  }
  if (aa == 't' || aa == 'T')
  {
    return 16;
  }
  if (aa == 'v' || aa == 'V')
  {
    return 17;
  }
  if (aa == 'w' || aa == 'W')
  {
    return 18;
  }
  if (aa == 'y' || aa == 'Y')
  {
    return 19;
  }
  return -1;
}

/* Converts a numerical value to an amino acid letter */
char amino_letter(int num)
{
  if (num == 0)
  {
    return 'A';
  }
  if (num == 1)
  {
    return 'C';
  }
  if (num == 2)
  {
    return 'D';
  }
  if (num == 3)
  {
    return 'E';
  }
  if (num == 4)
  {
    return 'F';
  }
  if (num == 5)
  {
    return 'G';
  }
  if (num == 6)
  {
    return 'H';
  }
  if (num == 7)
  {
    return 'I';
  }
  if (num == 8)
  {
    return 'K';
  }
  if (num == 9)
  {
    return 'L';
  }
  if (num == 10)
  {
    return 'M';
  }
  if (num == 11)
  {
    return 'N';
  }
  if (num == 12)
  {
    return 'P';
  }
  if (num == 13)
  {
    return 'Q';
  }
  if (num == 14)
  {
    return 'R';
  }
  if (num == 15)
  {
    return 'S';
  }
  if (num == 16)
  {
    return 'T';
  }
  if (num == 17)
  {
    return 'V';
  }
  if (num == 18)
  {
    return 'W';
  }
  if (num == 19)
  {
    return 'Y';
  }
  return 'X';
}

/* Assign the appropriate start value to this node */
int assign_start_value(unsigned char *seq, int n)
{
  if (is_atg(seq, n) == 1)
  {
    return ATG;
  }
  else if (is_gtg(seq, n) == 1)
  {
    return GTG;
  }
  else if (is_ttg(seq, n) == 1)
  {
    return TTG;
  }
  else
  {
    return NONST;
  }
}

/* Assign the appropriate stop value to this node */
int assign_stop_value(unsigned char *seq, int n)
{
  if (is_taa(seq, n) == 1)
  {
    return TAA;
  }
  else if (is_tag(seq, n) == 1)
  {
    return TAG;
  }
  else if (is_tga(seq, n) == 1)
  {
    return TGA;
  }
  else
  {
    return NONST;
  }
}

/* Returns the corresponding frame on the reverse strand */
int rev_frame(int frame, int seq_len)
{
  int mod = seq_len%3-1;
  if (mod == 0)
  {
    mod = 3;
  }
  return (mod-frame);
}

/* Simple 3-way max function for frames */
int max_frame(int n1, int n2, int n3)
{
  if (n1 > n2)
  {
    if (n1 > n3)
    {
      return 0;
    }
    else
    {
      return 2;
    }
  }
  else
  {
    if (n2 > n3)
    {
      return 1;
    }
    else
    {
      return 2;
    }
  }
}

/******************************************************************************
  Creates a GC frame plot for a given sequence.  This is simply a string with
  the highest GC content frame in a window centered on a position for every
  position in the sequence.
******************************************************************************/
int *calc_most_gc_frame(unsigned char *seq, int seq_length)
{
  int i = 0;
  int j = 0;
  int *for_sum = NULL;     /* Running gc_count of GC frame plot left->right */
  int *back_sum = NULL;    /* Running gc_count of GC frame plot right->left */
  int *gc_count = NULL;    /* GC count in window centered on position */
  int *gc_frame = NULL;    /* Frame with most GC in this position */
  int win = 0;

  gc_frame = (int *)calloc(seq_length, sizeof(int));
  for_sum = (int *)calloc(seq_length, sizeof(int));
  back_sum = (int *)calloc(seq_length, sizeof(int));
  gc_count = (int *)calloc(seq_length, sizeof(int));
  if (for_sum == NULL || back_sum == NULL || gc_frame == NULL ||
      gc_count == NULL)
  {
    perror("\n\nCalloc failed on GC frame plot.");
    exit(11);
  }
  for (i = 0; i < seq_length; i++)
  {
    gc_frame[i] = -1;
  }

  /* Calculate running gc_count GC counts */
  for (i = 0; i < 3; i++)
  {
    for (j = 0 + i; j < seq_length; j++)
    {
      if (j < 3)
      {
        for_sum[j] = is_gc(seq, j);
      }
      else
      {
        for_sum[j] = for_sum[j-3] + is_gc(seq, j);
      }
      if (j < 3)
      {
        back_sum[seq_length-j-1] = is_gc(seq, seq_length-j-1);
      }
      else
      {
        back_sum[seq_length-j-1] = back_sum[seq_length-j+2] +
                                   is_gc(seq, seq_length-j-1);
      }
    }
  }
  /* Using the counts, find the GC in a particular window */
  for (i = 0; i < seq_length; i++)
  {
    gc_count[i] = for_sum[i] + back_sum[i] - is_gc(seq, i);
    if (i - WINDOW/2 >= 0)
    {
      gc_count[i] -= for_sum[i-WINDOW/2];
    }
    if (i + WINDOW/2 < seq_length)
    {
      gc_count[i] -= back_sum[i+WINDOW/2];
    }
  }
  free(for_sum);
  free(back_sum);
  /* Find frame in window with highest GC count */
  for (i = 0; i < seq_length-2; i+=3)
  {
    win = max_frame(gc_count[i], gc_count[i+1], gc_count[i+2]);
    for (j = 0; j < 3; j++)
    {
      gc_frame[i+j] = win;
    }
  }
  free(gc_count);
  return gc_frame;
}

/* Converts a word of size len to a number 1 to 4^len */
int mer_index(int len, unsigned char *seq, int pos)
{
  int i = 0;
  int index = 0;
  for (i = 0; i < 2*len; i++)
  {
    index |= (test(seq, pos*2+i)<<i);
  }
  return index;
}

/* Gives a text string for a start */
void start_text(char *st, int type)
{
  if (type == 0)
  {
    st[0] = 'A';
  }
  else if (type == 1)
  {
    st[0] = 'G';
  }
  else if (type == 2)
  {
    st[0] = 'T';
  }
  st[1] = 'T';
  st[2] = 'G';
  st[3] = '\0';
}

/* Gives a text string for a mer of len 'len' (useful in outputting motifs) */
void mer_text(char *text, int len, int index)
{
  int i = 0;
  int val = 0;      /* Numerical value of the letter */
  char letters[4] = { 'A', 'G', 'C', 'T' };
  if (len == 0)
  {
    strcpy(text, "None");
  }
  else
  {
    for (i = 0; i < len; i++)
    {
      val = (index&(1<<(2*i))) + (index&(1<<(2*i+1)));
      val >>= (i*2);
      text[i] = letters[val];
    }
    text[i] = '\0';
  }
}

/* Builds a 'len'-mer background for whole sequence */
void calc_mer_background(int len, unsigned char *seq, unsigned char *rseq,
                         int seq_length, double *background)
{
  int i = 0;
  int sum = 0;          /* Sum of all counts */
  int size = 1;         /* Size is 4^len */
  int *counts = NULL;   /* Counts of particular words */

  for (i = 1; i <= len; i++)
  {
    size *= 4;
  }
  counts = (int *)calloc(size, sizeof(int));
  if (counts == NULL)
  {
    perror("\n\nCalloc failed on mer background.");
    exit(11);
  }
  for (i = 0; i < seq_length-len+1; i++)
  {
    counts[mer_index(len, seq, i)]++;
    counts[mer_index(len, rseq, i)]++;
    sum += 2;
  }
  for (i = 0; i < size; i++)
  {
    background[i] = (double)((counts[i]*1.0)/(sum*1.0));
  }
  free(counts);
}

/******************************************************************************
  For a given start, record the base composition of the -1/-2 and -15 to -44
  positions upstream.
******************************************************************************/
void count_upstream_composition(unsigned char *seq, int seq_length, int strand,
                                int pos, int *ups)
{
  int i = 0;
  int start = 0;              /* Start site position */
  int counter = 0;            /* Counter */

  if (strand == 1)
  {
    start = pos;
  }
  else
  {
    start = seq_length-1-pos;
  }
  if (start - 45 < 0)
  {
    return;
  }

  for (i = start-44; i < start; i++)
  {
    if (i >= start-14 && i <= start-3)
    {
      continue;
    }
    ups[counter] = mer_index(1, seq, i);
    counter++;
  }
}

/******************************************************************************
  Finds the highest-scoring subset of AGGAGG in a given stretch of
  sequence upstream of a start (the sequence between ups_start and
  gene_start), using rbs_weight values to score the motifs in each bin.
******************************************************************************/
int shine_dalgarno_exact(unsigned char *seq, int ups_start, int gene_start,
                         double *rbs_weight)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int mismatch = 0;             /* To track if we see a mismatch */
  int spacer = 0;               /* Distance between end of motif and gene */
  int spacer_bin = 0;           /* Convert spacer to various bins */
  int limit = 0;                /* Used to handle really close ups_start */
  int max_bin = 0;              /* Maximum motif bin we've seen so far */
  int motif_bin = 0;            /* Current motif bin */
  double match[6] = {0};        /* Scoring (used to assign bins) */
  double match_sum = 0.0;       /* Total score of this motif */

  limit = imin(6, gene_start-4-ups_start);
  for (i = limit; i < 6; i++)
  {
    match[i] = -10.0;
  }

  /* Compare the 6-base region to AGGAGG */
  for (i = 0; i < limit; i++)
  {
    if (ups_start+i < 0)
    {
      continue;
    }
    if (i%3 == 0 && is_a(seq, ups_start+i) == 1)
    {
      match[i] = 2.0;
    }
    else if (i%3 != 0 && is_g(seq, ups_start+i) == 1)
    {
      match[i] = 3.0;
    }
    else
    {
      match[i] = -10.0;
    }
  }

  /* Find the maximally scoring motif */
  max_bin = 0;
  for (i = limit; i >= 3; i--)
  {
    for (j = 0; j <= limit-i; j++)
    {
      match_sum = -2.0;
      mismatch = 0;
      for (k = j; k < j+i; k++)
      {
        match_sum += match[k];
        if (match[k] < 0.0)
        {
          mismatch++;
        }
      }
      if (mismatch > 0)
      {
        continue;
      }
      /* Calculate spacer and assign it to a bin */
      spacer = gene_start - (ups_start+j+i);
      if (spacer < 5 && i < 5)
      {
        spacer_bin = 2;
      }
      else if (spacer < 5 && i >= 5)
      {
        spacer_bin = 1;
      }
      else if (spacer > 10 && spacer <= 12 && i < 5)
      {
        spacer_bin = 1;
      }
      else if (spacer > 10 && spacer <= 12 && i >= 5)
      {
        spacer_bin = 2;
      }
      else if (spacer >= 13)
      {
        spacer_bin = 3;
      }
      else
      {
        spacer_bin = 0;
      }
      if (spacer > 15 || match_sum < 6.0)
      {
        continue;
      }

      /* Exact-Matching RBS Motifs */
      /* Map the match score and spacer bin to an overall motif bin */
      if (match_sum < 6.0)
      {
        motif_bin = 0;
      }
      else if (match_sum == 6.0 && spacer_bin == 2)
      {
        motif_bin = 1;
      }
      else if (match_sum == 6.0 && spacer_bin == 3)
      {
        motif_bin = 2;
      }
      else if (match_sum == 8.0 && spacer_bin == 3)
      {
        motif_bin = 3;
      }
      else if (match_sum == 9.0 && spacer_bin == 3)
      {
        motif_bin = 3;
      }
      else if (match_sum == 6.0 && spacer_bin == 1)
      {
        motif_bin = 6;
      }
      else if (match_sum == 11.0 && spacer_bin == 3)
      {
        motif_bin = 10;
      }
      else if (match_sum == 12.0 && spacer_bin == 3)
      {
        motif_bin = 10;
      }
      else if (match_sum == 14.0 && spacer_bin == 3)
      {
        motif_bin = 10;
      }
      else if (match_sum == 8.0 && spacer_bin == 2)
      {
        motif_bin = 11;
      }
      else if (match_sum == 9.0 && spacer_bin == 2)
      {
        motif_bin = 11;
      }
      else if (match_sum == 8.0 && spacer_bin == 1)
      {
        motif_bin = 12;
      }
      else if (match_sum == 9.0 && spacer_bin == 1)
      {
        motif_bin = 12;
      }
      else if (match_sum == 6.0 && spacer_bin == 0)
      {
        motif_bin = 13;
      }
      else if (match_sum == 8.0 && spacer_bin == 0)
      {
        motif_bin = 15;
      }
      else if (match_sum == 9.0 && spacer_bin == 0)
      {
        motif_bin = 16;
      }
      else if (match_sum == 11.0 && spacer_bin == 2)
      {
        motif_bin = 20;
      }
      else if (match_sum == 11.0 && spacer_bin == 1)
      {
        motif_bin = 21;
      }
      else if (match_sum == 11.0 && spacer_bin == 0)
      {
        motif_bin = 22;
      }
      else if (match_sum == 12.0 && spacer_bin == 2)
      {
        motif_bin = 20;
      }
      else if (match_sum == 12.0 && spacer_bin == 1)
      {
        motif_bin = 23;
      }
      else if (match_sum == 12.0 && spacer_bin == 0)
      {
        motif_bin = 24;
      }
      else if (match_sum == 14.0 && spacer_bin == 2)
      {
        motif_bin = 25;
      }
      else if (match_sum == 14.0 && spacer_bin == 1)
      {
        motif_bin = 26;
      }
      else if (match_sum == 14.0 && spacer_bin == 0)
      {
        motif_bin = 27;
      }

      /* Use the values passed in by rbs_weight to find the score */
      /* for this bin.  For initial runs, when rbs_weights are all 0, */
      /* the highest numbered bin wins. */
      if (rbs_weight[motif_bin] > rbs_weight[max_bin] ||
          (rbs_weight[motif_bin] == rbs_weight[max_bin] &&
          motif_bin >= max_bin))
      {
        max_bin = motif_bin;
      }
    }
  }
  return max_bin;
}

/******************************************************************************
  Finds the highest-scoring subset of AGGAGG with at least 1 mismatch in a
  given stretch of sequence upstream of a start (the sequence between
  ups_start and gene_start), using rbs_weight values to score the motifs in
  each bin.
******************************************************************************/
int shine_dalgarno_mismatch(unsigned char *seq, int ups_start, int gene_start,
                            double *rbs_weight)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int mismatch = 0;             /* To track if we see a mismatch */
  int spacer = 0;               /* Distance between end of motif and gene */
  int spacer_bin = 0;           /* Convert spacer to various bins */
  int limit = 0;                /* Used to handle really close ups_start */
  int max_bin = 0;              /* Maximum motif bin we've seen so far */
  int motif_bin = 0;            /* Current motif bin */
  double match[6] = {0};        /* Scoring (used to assign bins) */
  double match_sum = 0.0;       /* Total score of this motif */

  limit = imin(6, gene_start-4-ups_start);
  for (i = limit; i < 6; i++)
  {
    match[i] = -10.0;
  }

  /* Compare the 6-base region to AGGAGG */
  for (i = 0; i < limit; i++)
  {
    if (ups_start+i < 0)
    {
      continue;
    }
    if (i % 3 == 0)
    {
      if (is_a(seq, ups_start+i) == 1)
      {
        match[i] = 2.0;
      }
      else
      {
        match[i] = -3.0;
      }
    }
    else
    {
      if (is_g(seq, ups_start+i) == 1)
      {
        match[i] = 3.0;
      }
      else
      {
        match[i] = -2.0;
      }
    }
  }

  /* Find the maximally scoring motif */
  max_bin = 0;
  for (i = limit; i >= 5; i--)
  {
    for (j = 0; j <= limit-i; j++)
    {
      match_sum = -2.0;
      mismatch = 0;
      for (k = j; k < j+i; k++)
      {
        match_sum += match[k];
        if (match[k] < 0.0)
        {
          mismatch++;
        }
        if (match[k] < 0.0 && (k <= j+1 || k >= j+i-2))
        {
          match_sum -= 10.0;
        }
      }
      if (mismatch != 1)
      {
        continue;
      }
      spacer = gene_start - (ups_start+j+i);
      if (spacer < 5)
      {
        spacer_bin = 1;
      }
      else if (spacer > 10 && spacer <= 12)
      {
        spacer_bin = 2;
      }
      else if (spacer >= 13)
      {
        spacer_bin = 3;
      }
      else
      {
        spacer_bin = 0;
      }
      if (spacer > 15 || match_sum < 6.0)
      {
        continue;
      }

      /* Single-Mismatch RBS Motifs */
      if (match_sum < 6.0)
      {
        motif_bin = 0;
      }
      else if (match_sum == 6.0 && spacer_bin == 3)
      {
        motif_bin = 2;
      }
      else if (match_sum == 7.0 && spacer_bin == 3)
      {
        motif_bin = 2;
      }
      else if (match_sum == 9.0 && spacer_bin == 3)
      {
        motif_bin = 3;
      }
      else if (match_sum == 6.0 && spacer_bin == 2)
      {
        motif_bin = 4;
      }
      else if (match_sum == 6.0 && spacer_bin == 1)
      {
        motif_bin = 5;
      }
      else if (match_sum == 6.0 && spacer_bin == 0)
      {
        motif_bin = 9;
      }
      else if (match_sum == 7.0 && spacer_bin == 2)
      {
        motif_bin = 7;
      }
      else if (match_sum == 7.0 && spacer_bin == 1)
      {
        motif_bin = 8;
      }
      else if (match_sum == 7.0 && spacer_bin == 0)
      {
        motif_bin = 14;
      }
      else if (match_sum == 9.0 && spacer_bin == 2)
      {
        motif_bin = 17;
      }
      else if (match_sum == 9.0 && spacer_bin == 1)
      {
        motif_bin = 18;
      }
      else if (match_sum == 9.0 && spacer_bin == 0)
      {
        motif_bin = 19;
      }

      if (rbs_weight[motif_bin] < rbs_weight[max_bin])
      {
        continue;
      }
      if (rbs_weight[motif_bin] == rbs_weight[max_bin] && motif_bin < max_bin)
      {
        continue;
      }
      max_bin = motif_bin;
    }
  }

  return max_bin;
}

/* Zero out the sequence bitstrings for reuse */
void zero_sequence(unsigned char *seq, unsigned char *rseq,
                   unsigned char *useq, int slen)
{
  memset(seq, 0, (slen/4 + 1) * sizeof(unsigned char));
  memset(rseq, 0, (slen/4 + 1) * sizeof(unsigned char));
  memset(useq, 0, (slen/8 + 1) * sizeof(unsigned char));
}

/* Returns the minimum of two integers */
int imin(int x, int y)
{
  if (x < y)
  {
    return x;
  }
  return y;
}

/* Returns the maximum of two integers */
int imax(int x, int y)
{
  if (x > y)
  {
    return x;
  }
  return y;
}

/* Return the minimum of two numbers */
double dmin(double x, double y)
{
  if (x < y)
  {
    return x;
  }
  return y;
}

/* Return the maximum of two numbers */
double dmax(double x, double y)
{
  if (x > y)
  {
    return x;
  }
  return y;
}
