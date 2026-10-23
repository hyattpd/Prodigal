/*******************************************************************************
    File and FASTA input for the test programs (test_api.c and tests/),
    static inline so each uses what it needs. They abort when memory runs
    out.

    fasta_load keeps the records in order, so record k is the CLI's sequence
    number k+1. Line ends (LF or CRLF) are stripped; every other byte is
    passed through, because prodigal_set_sequence applies the same character
    rules as the CLI's reader. A record with no sequence lines has length 0.
*******************************************************************************/

#ifndef PRODIGAL_TESTS_FASTA_H
#define PRODIGAL_TESTS_FASTA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static inline void *xrealloc(void *p, size_t n) {
    if ((p = realloc(p, n ? n : 1)) == NULL) {
        fputs("out of memory\n", stderr);
        abort();
    }
    return p;
}

static inline char *dupstr(const char *s) {
    size_t n = strlen(s) + 1;
    return memcpy(xrealloc(NULL, n), s, n);
}

/* The contents of the file at path, or NULL if it can't be read. */
static inline void *load_file(const char *path, size_t *len) {
    FILE *fp = fopen(path, "rb");
    char *buf = NULL;
    long n;
    if (fp == NULL) return NULL;
    if (fseek(fp, 0, SEEK_END) == 0 && (n = ftell(fp)) >= 0 && fseek(fp, 0, SEEK_SET) == 0) {
        buf = xrealloc(NULL, (size_t)n);
        if (fread(buf, 1, (size_t)n, fp) == (size_t)n) *len = (size_t)n;
        else { free(buf); buf = NULL; }
    }
    fclose(fp);
    return buf;
}

typedef struct {
    int32_t n;
    char **seqs;     /* n sequences, each NUL-terminated */
    char **headers;  /* n headers without the '>' */
    int32_t *lens;   /* n lengths */
} fasta_t;

static inline void fasta_free(fasta_t *f) {
    int32_t i;
    for (i = 0; i < f->n; i++) { free(f->seqs[i]); free(f->headers[i]); }
    free(f->seqs); free(f->headers); free(f->lens);
    memset(f, 0, sizeof(*f));
}

#define FASTA_LINE_MAX 65536

/* Returns 0, or -1 if the file can't be opened. */
static inline int fasta_load(const char *path, fasta_t *f) {
    FILE *fp = fopen(path, "rb");
    char *line;
    int32_t cap = 0, scap = 0;

    memset(f, 0, sizeof(*f));
    if (fp == NULL) return -1;
    line = xrealloc(NULL, FASTA_LINE_MAX);   /* not on the stack: WASM's is small */
    while (fgets(line, FASTA_LINE_MAX, fp) != NULL) {
        size_t ln = strlen(line);
        while (ln > 0 && (line[ln - 1] == '\n' || line[ln - 1] == '\r'))
            line[--ln] = '\0';
        if (line[0] == '>') {
            int32_t k = f->n++;
            if (k == cap) {
                cap = cap ? cap * 2 : 64;
                f->seqs = xrealloc(f->seqs, cap * sizeof(char *));
                f->headers = xrealloc(f->headers, cap * sizeof(char *));
                f->lens = xrealloc(f->lens, cap * sizeof(int32_t));
            }
            f->headers[k] = dupstr(line + 1);
            f->seqs[k] = dupstr("");
            f->lens[k] = 0;
            scap = 1;
        } else if (f->n > 0 && ln > 0) {
            int32_t k = f->n - 1, need = f->lens[k] + (int32_t)ln + 1;
            if (need > scap) {
                while (scap < need) scap *= 2;
                f->seqs[k] = xrealloc(f->seqs[k], scap);
            }
            memcpy(f->seqs[k] + f->lens[k], line, ln + 1);
            f->lens[k] += (int32_t)ln;
        }
    }
    fclose(fp);
    free(line);
    return 0;
}

#endif
