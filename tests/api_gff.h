/*******************************************************************************
    CLI-format lines from the public API's output, for comparing the library
    with the CLI (test_api.c, tests/api_dump.c, tests/threads.c); static
    inline, so each uses what it needs.

    Per contig, a model line, then one line per gene:
      #  seqnum  seqlen=N  model=DESC  transl_table=T  uses_sd=U
      seqid  CDS  begin  end  score  strand  0  ID=..;partial=..;...;tscore=..;
    The gene line is the CLI's GFF CDS line without column 2 (the version);
    the model line holds the fields of its "# Sequence Data" / "# Model Data"
    lines that the API reports (model is "-" in single mode). Numbers use the
    CLI's formats, so equal values print identically. When the genes carry
    proteins (want_proteins), each gene line ends with a tab and the protein;
    on the CLI side, the records of its -a file, in order.

    Column 6 is "." when blank_score is set. In metagenomic mode the CLI
    prints it from its recovery pass rather than from the values that chose
    the gene (docs/upstream-behavior.md, quirk 2), so it is not part of
    the API contract; the score= attribute is.
*******************************************************************************/

#ifndef PRODIGAL_TESTS_API_GFF_H
#define PRODIGAL_TESTS_API_GFF_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../prodigal.h"
#include "fasta.h"

#define API_GFF_LINE_MAX 2048

typedef void (*api_emit_fn)(const char *line, void *arg);

/* Check snprintf's result r for a buffer of n bytes: a line that doesn't fit
   is a harness bug, never a silent difference. */
static inline void api_fits(int r, size_t n) {
    if (r < 0 || (size_t)r >= n) {
        fputs("api_gff: line too long\n", stderr);
        abort();
    }
}

/* The CLI's sequence ID: the header up to the first space, tab or CR, or
   "Prodigal_Seq_<seqnum>" when that is empty (calc_short_header). */
static inline void api_seqid(char *out, size_t n, const char *hdr, int seqnum) {
    size_t i = strcspn(hdr, " \t\r\n");
    if (i == 0) snprintf(out, n, "Prodigal_Seq_%d", seqnum);
    else snprintf(out, n, "%.*s", (int)i, hdr);
}

/* The model line. model ("-" in single mode) ends at a NUL or a double
   quote, so the CLI side can pass the rest of its line. */
static inline void api_model_line(char *out, size_t n, int seqnum, int seqlen,
                                  const char *model, int table, int uses_sd) {
    api_fits(snprintf(out, n, "#\t%d\tseqlen=%d\tmodel=%.*s\ttransl_table=%d\tuses_sd=%d",
                      seqnum, seqlen, (int)strcspn(model, "\""), model, table,
                      uses_sd), n);
}

/* Gene i of g. seqid is api_seqid's, seqnum the 1-based record index. */
static inline void api_gff_line(char *out, size_t n, const prodigal_genes_soa_t *g,
                                int i, const char *seqid, int seqnum,
                                int blank_score) {
    static const char *const start_types[4] = { "ATG", "GTG", "TTG", "Edge" };
    double score = g->cscore[i] + g->sscore[i];
    char col6[32];
    int st = g->start_type[i];

    if (blank_score) strcpy(col6, ".");
    else snprintf(col6, sizeof(col6), "%.1f", score);
    api_fits(snprintf(out, n,
                 "%s\tCDS\t%d\t%d\t%s\t%c\t0\t"
                 "ID=%d_%d;partial=%d%d;start_type=%s;"
                 "rbs_motif=%s;rbs_spacer=%s;gc_cont=%.3f;"
                 "conf=%.2f;score=%.2f;cscore=%.2f;sscore=%.2f;"
                 "rscore=%.2f;uscore=%.2f;tscore=%.2f;",
                 seqid, g->begin[i], g->end[i], col6,
                 g->strand[i] == 1 ? '+' : '-',
                 seqnum, i + 1, g->partial_left[i], g->partial_right[i],
                 st >= 0 && st < 4 ? start_types[st] : "?",
                 g->rbs_motif[i], g->rbs_spacer[i], g->gc_cont[i],
                 g->confidence[i], score, g->cscore[i], g->sscore[i],
                 g->rscore[i], g->uscore[i], g->tscore[i]), n);
}

/* Gene i's protein in g (want_proteins): its letters and, in len, how
   many. */
static inline const char *api_protein(const prodigal_genes_soa_t *g, int32_t i,
                                      int32_t *len) {
    *len = g->protein_offsets[i + 1] - g->protein_offsets[i];
    return g->protein_data + g->protein_offsets[i];
}

/* 1 if g's proteins are well formed: aligned offsets from 0, a letter per
   codon, M at a start codon (start type not Edge), and * only at a stop
   codon, none when the gene runs off the edge at its stop end (where a
   table-22 TNA stop gives X: docs/upstream-behavior.md, quirk 13).
   Otherwise 0, with what is wrong in why. */
static inline int api_proteins_ok(const prodigal_genes_soa_t *g, char *why,
                                  size_t n) {
    int32_t i, k, len;
    if (g->n_genes == 0) return 1;
    if (g->protein_offsets == NULL || g->protein_data == NULL ||
        ((uintptr_t)g->protein_offsets) % 16 != 0 || g->protein_offsets[0] != 0) {
        snprintf(why, n, "protein arrays missing, misaligned or not from 0");
        return 0;
    }
    for (i = 0; i < g->n_genes; i++) {
        const char *p = api_protein(g, i, &len);
        int stop_partial = g->strand[i] == 1 ? g->partial_right[i] : g->partial_left[i];
        if (len != (g->end[i] - g->begin[i] + 1) / 3) {
            snprintf(why, n, "gene %d: %d letters for %d..%d", i + 1, len,
                     g->begin[i], g->end[i]);
            return 0;
        }
        for (k = 0; k < len; k++) {
            int stop = k == len - 1 && !stop_partial;
            if (p[k] == '\0' || strchr("ACDEFGHIKLMNPQRSTVWXY*", p[k]) == NULL ||
                (stop ? p[k] != '*' && p[k] != 'X' : p[k] == '*') ||
                (k == 0 && g->start_type[i] != 3 && p[k] != 'M')) {
                snprintf(why, n, "gene %d: letter %d of %d is '%c'", i + 1, k + 1,
                         len, p[k] ? p[k] : '0');
                return 0;
            }
        }
    }
    return 1;
}

/* Emit line, then a tab and the n letters of prot. */
static inline void api_emit_with_protein(const char *line, const char *prot,
                                         size_t n, api_emit_fn emit, void *arg) {
    size_t ln = strlen(line);
    char *buf = (char *)xrealloc(NULL, ln + n + 2);
    memcpy(buf, line, ln);
    buf[ln] = '\t';
    memcpy(buf + ln + 1, prot, n);
    buf[ln + 1 + n] = '\0';
    emit(buf, arg);
    free(buf);
}

/* Reader state carried from a "# Sequence Data" line to its "# Model Data". */
typedef struct { int seqnum, seqlen; } api_cli_state_t;

/* Value of key in line ("key=VALUE;..."), or NULL. */
static inline const char *api_cli_field(const char *line, const char *key) {
    const char *p = strstr(line, key);
    return p ? p + strlen(key) : NULL;
}

/* Translate one line of a CLI GFF file into the form above; line is split in
   place. Fills out and returns API_CLI_MODEL or API_CLI_GENE for a model or
   CDS line; returns -1 for anything else. */
enum { API_CLI_GENE, API_CLI_MODEL };
static inline int api_gff_from_cli(char *line, int blank_score, api_cli_state_t *st,
                                   char *out, size_t n) {
    char *col[9];
    int k;

    if (strncmp(line, "# Sequence Data: ", 17) == 0) {
        const char *a = api_cli_field(line, "seqnum="), *b = api_cli_field(line, ";seqlen=");
        st->seqnum = a ? atoi(a) : -1;
        st->seqlen = b ? atoi(b) : -1;
        return -1;
    }
    if (strncmp(line, "# Model Data: ", 14) == 0) {
        const char *rt = api_cli_field(line, "run_type=");
        const char *m = api_cli_field(line, "model=\"");
        const char *tt = api_cli_field(line, "transl_table=");
        const char *sd = api_cli_field(line, "uses_sd=");
        if (rt == NULL || strncmp(rt, "Metagenomic", 11) != 0 || m == NULL) m = "-";
        api_model_line(out, n, st->seqnum, st->seqlen, m, tt ? atoi(tt) : -1,
                       sd ? atoi(sd) : -1);
        return API_CLI_MODEL;
    }
    if (line[0] == '#') return -1;

    line[strcspn(line, "\r\n")] = '\0';
    for (col[0] = line, k = 1; k < 9; k++) {
        if ((col[k] = strchr(col[k - 1], '\t')) == NULL) return -1;
        *col[k]++ = '\0';
    }
    if (strcmp(col[2], "CDS") != 0) return -1;
    api_fits(snprintf(out, n, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", col[0], col[2],
                      col[3], col[4], blank_score ? "." : col[5], col[6], col[7],
                      col[8]), n);
    return API_CLI_GENE;
}

/* Train ctx on all of fa, as the CLI does in single mode. Returns
   PRODIGAL_OK or the failing call's error code. */
static inline int api_train(prodigal_ctx_t *ctx, const fasta_t *fa) {
    int rc = prodigal_set_training_sequences(ctx, (const char **)fa->seqs,
                 (const char **)fa->headers, fa->lens, fa->n);
    if (rc == PRODIGAL_OK) rc = prodigal_train(ctx);
    return rc;
}

/* Find genes in the sequence set on ctx, record k of fa, and emit its model
   line and gene lines; meta selects metagenomic formatting. Malformed
   proteins (api_proteins_ok) add a line "BAD proteins: <why>", so a
   comparison fails on them. Returns PRODIGAL_OK or prodigal_find_genes's
   error code. */
static inline int api_emit_found(prodigal_ctx_t *ctx, const fasta_t *fa, int32_t k,
                                 int meta, api_emit_fn emit, void *arg) {
    char seqid[API_GFF_LINE_MAX], line[API_GFF_LINE_MAX], why[128];
    prodigal_genes_soa_t *g = NULL;
    prodigal_stats_t stats;
    prodigal_seq_info_t info;
    int32_t j, len;
    int rc = prodigal_find_genes(ctx, &g, &stats);
    if (rc != PRODIGAL_OK) return rc;
    prodigal_get_seq_info(ctx, &info);
    api_model_line(line, sizeof(line), k + 1, info.length,
                   meta ? stats.best_meta_desc : "-",
                   stats.translation_table, stats.uses_sd);
    emit(line, arg);
    api_seqid(seqid, sizeof(seqid), fa->headers[k], k + 1);
    for (j = 0; j < g->n_genes; j++) {
        api_gff_line(line, sizeof(line), g, j, seqid, k + 1, meta);
        if (g->protein_offsets == NULL) {
            emit(line, arg);
        } else {
            const char *p = api_protein(g, j, &len);
            api_emit_with_protein(line, p, (size_t)len, emit, arg);
        }
    }
    if (g->protein_offsets != NULL && !api_proteins_ok(g, why, sizeof(why))) {
        snprintf(line, sizeof(line), "BAD proteins: %s", why);
        emit(line, arg);
    }
    prodigal_genes_free(g);
    return PRODIGAL_OK;
}

/* api_emit_found for records [first, last) of fa, each set on ctx first.
   Returns PRODIGAL_OK or the failing call's error code. */
static inline int api_emit_records(prodigal_ctx_t *ctx, const fasta_t *fa,
                                   int32_t first, int32_t last, int meta,
                                   api_emit_fn emit, void *arg) {
    int32_t k;
    for (k = first; k < last; k++) {
        int rc = prodigal_set_sequence(ctx, fa->seqs[k], fa->lens[k], fa->headers[k]);
        if (rc == PRODIGAL_OK) rc = api_emit_found(ctx, fa, k, meta, emit, arg);
        if (rc != PRODIGAL_OK) return rc;
    }
    return PRODIGAL_OK;
}

/* Emitted lines, collected (lines_emit) and compared (lines_equal). */
typedef struct { char **v; int n, cap; } lines_t;

static inline void lines_emit(const char *line, void *arg) {
    lines_t *l = (lines_t *)arg;
    if (l->n == l->cap) {
        l->cap = l->cap ? l->cap * 2 : 256;
        l->v = (char **)xrealloc(l->v, l->cap * sizeof(char *));
    }
    l->v[l->n++] = dupstr(line);
}

static inline void lines_free(lines_t *l) {
    int i;
    for (i = 0; i < l->n; i++) free(l->v[i]);
    free(l->v);
    memset(l, 0, sizeof(*l));
}

/* 1 if want and got are equal; otherwise prints [FAIL] with the line counts,
   the number of differing lines, and the first difference (or extra line),
   and returns 0. */
static inline int lines_equal(const lines_t *want, const lines_t *got,
                              const char *what) {
    int n = want->n < got->n ? want->n : got->n, i, ndiff = 0, first = -1;
    for (i = 0; i < n; i++)
        if (strcmp(want->v[i], got->v[i]) != 0 && ndiff++ == 0) first = i;
    if (ndiff == 0 && want->n == got->n) return 1;
    if (first < 0) first = n;
    printf("[FAIL] %s: %d lines (expected %d), %d of the first %d differ; "
           "first, line %d:\n      want %s\n      got  %s\n", what, got->n,
           want->n, ndiff, n, first + 1,
           first < want->n ? want->v[first] : "(none)",
           first < got->n ? got->v[first] : "(none)");
    return 0;
}

/* Emit a CLI GFF file in the same form, with the proteins of its -a file
   faa, in order, unless that is NULL. Returns 0, API_CLI_UNREADABLE if a
   file can't be read or the GFF has a line too long for the buffer, or
   API_CLI_MISMATCH if the two files hold different numbers of genes. */
enum { API_CLI_UNREADABLE = -1, API_CLI_MISMATCH = -2 };
static inline int api_emit_cli_file(const char *path, const char *faa, int meta,
                                    api_emit_fn emit, void *arg) {
    char line[API_GFF_LINE_MAX], out[API_GFF_LINE_MAX];
    api_cli_state_t st = { -1, -1 };
    fasta_t prot = {0};
    FILE *fp;
    int32_t ngenes = 0;
    int rc = 0, kind;
    if (faa != NULL && fasta_load(faa, &prot) != 0) return API_CLI_UNREADABLE;
    if ((fp = fopen(path, "r")) == NULL) {
        fasta_free(&prot);
        return API_CLI_UNREADABLE;
    }
    while (rc == 0 && fgets(line, sizeof(line), fp) != NULL) {
        if (strchr(line, '\n') == NULL && !feof(fp)) rc = API_CLI_UNREADABLE;
        else if ((kind = api_gff_from_cli(line, meta, &st, out, sizeof(out))) < 0) continue;
        else if (faa == NULL || kind == API_CLI_MODEL) emit(out, arg);
        else if (ngenes == prot.n) rc = API_CLI_MISMATCH;
        else {
            api_emit_with_protein(out, prot.seqs[ngenes], (size_t)prot.lens[ngenes],
                                  emit, arg);
            ngenes++;
        }
    }
    if (rc == 0)
        rc = !feof(fp) ? API_CLI_UNREADABLE : ngenes != prot.n ? API_CLI_MISMATCH : 0;
    fclose(fp);
    fasta_free(&prot);
    return rc;
}

#endif
