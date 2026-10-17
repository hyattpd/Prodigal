/*******************************************************************************
    CLI-format lines from the public API's output, for comparing the library
    with the CLI (test_api.c, tests/api_dump.c).

    Per contig, a model line, then one line per gene:
      #  seqnum  seqlen=N  model=DESC  transl_table=T  uses_sd=U
      seqid  CDS  begin  end  score  strand  0  ID=..;partial=..;...;tscore=..;
    The gene line is the CLI's GFF CDS line without column 2 (the version);
    the model line holds the fields of its "# Sequence Data" / "# Model Data"
    lines that the API reports (model is "-" in single mode). Numbers use the
    CLI's formats, so equal values print identically.

    Column 6 is "." when blank_score is set. In metagenomic mode the CLI
    prints it from its recovery pass rather than from the values that chose
    the gene (docs/upstream-behavior.md, quirk 2), so it is not part of
    the API contract; the score= attribute is.
*******************************************************************************/

#ifndef PRODIGAL_TESTS_API_GFF_H
#define PRODIGAL_TESTS_API_GFF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prodigal.h"
#include "fasta.h"

#define API_GFF_LINE_MAX 2048

/* Check snprintf's result r for a buffer of n bytes: a line that doesn't fit
   is a harness bug, never a silent difference. */
static void api_fits(int r, size_t n) {
    if (r < 0 || (size_t)r >= n) {
        fputs("api_gff: line too long\n", stderr);
        abort();
    }
}

/* The CLI's sequence ID: the header up to the first space, tab or CR, or
   "Prodigal_Seq_<seqnum>" when that is empty (calc_short_header). */
static void api_seqid(char *out, size_t n, const char *hdr, int seqnum) {
    size_t i = strcspn(hdr, " \t\r\n");
    if (i == 0) snprintf(out, n, "Prodigal_Seq_%d", seqnum);
    else snprintf(out, n, "%.*s", (int)i, hdr);
}

/* The model line. model ("-" in single mode) ends at a NUL or a double
   quote, so the CLI side can pass the rest of its line. */
static void api_model_line(char *out, size_t n, int seqnum, int seqlen,
                           const char *model, int table, int uses_sd) {
    api_fits(snprintf(out, n, "#\t%d\tseqlen=%d\tmodel=%.*s\ttransl_table=%d\tuses_sd=%d",
                      seqnum, seqlen, (int)strcspn(model, "\""), model, table,
                      uses_sd), n);
}

/* Gene i of g. seqid is api_seqid's, seqnum the 1-based record index. */
static void api_gff_line(char *out, size_t n, const prodigal_genes_soa_t *g,
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

/* Reader state carried from a "# Sequence Data" line to its "# Model Data". */
typedef struct { int seqnum, seqlen; } api_cli_state_t;

/* Value of key in line ("key=VALUE;..."), or NULL. */
static const char *api_cli_field(const char *line, const char *key) {
    const char *p = strstr(line, key);
    return p ? p + strlen(key) : NULL;
}

/* Replace the value of rbs_motif= in line, if any, with "?". */
static void api_mask_motif(char *line) {
    char *p = (char *)api_cli_field(line, "rbs_motif="), *e;
    if (p == NULL || (e = strchr(p, ';')) == NULL) return;
    *p = '?';
    memmove(p + 1, e, strlen(e) + 1);
}

/* Translate one line of a CLI GFF file into the form above; line is split in
   place. Returns 0 and fills out for a model or CDS line, -1 for anything
   else. */
static int api_gff_from_cli(char *line, int blank_score, api_cli_state_t *st,
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
        return 0;
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
    return 0;
}

/* Train ctx on all of fa, as the CLI does in single mode. fix_gc first sets
   the training GC by hand, which neutralises defect A (docs/testing.md).
   Returns PRODIGAL_OK or the failing call's error code. */
static int api_train(prodigal_ctx_t *ctx, const fasta_t *fa, int fix_gc) {
    prodigal_seq_info_t info;
    int rc = prodigal_set_training_sequences(ctx, (const char **)fa->seqs,
                 (const char **)fa->headers, fa->lens, fa->n);
    if (rc == PRODIGAL_OK && fix_gc) rc = prodigal_get_seq_info(ctx, &info);
    if (rc == PRODIGAL_OK && fix_gc) rc = prodigal_set_gc(ctx, info.gc_content);
    if (rc == PRODIGAL_OK) rc = prodigal_train(ctx);
    return rc;
}

typedef void (*api_emit_fn)(const char *line, void *arg);

/* Find genes in the sequence set on ctx, record k of fa, and emit its model
   line and gene lines; meta selects metagenomic formatting. Returns
   PRODIGAL_OK or prodigal_find_genes's error code. */
static int api_emit_found(prodigal_ctx_t *ctx, const fasta_t *fa, int32_t k,
                          int meta, api_emit_fn emit, void *arg) {
    char seqid[API_GFF_LINE_MAX], line[API_GFF_LINE_MAX];
    prodigal_genes_soa_t *g = NULL;
    prodigal_stats_t stats;
    prodigal_seq_info_t info;
    int32_t j;
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
        emit(line, arg);
    }
    prodigal_genes_free(g);
    return PRODIGAL_OK;
}

/* api_emit_found for records [first, last) of fa, each set on ctx first.
   Returns PRODIGAL_OK or the failing call's error code. */
static int api_emit_records(prodigal_ctx_t *ctx, const fasta_t *fa,
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

/* Emit a CLI GFF file in the same form; -1 if it can't be read or has a
   line too long for the buffer. */
static int api_emit_cli_file(const char *path, int meta, api_emit_fn emit,
                             void *arg) {
    char line[API_GFF_LINE_MAX], out[API_GFF_LINE_MAX];
    api_cli_state_t st = { -1, -1 };
    FILE *fp = fopen(path, "r");
    int k;
    if (fp == NULL) return -1;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (strchr(line, '\n') == NULL && !feof(fp)) break;
        if (api_gff_from_cli(line, meta, &st, out, sizeof(out)) == 0) emit(out, arg);
    }
    k = feof(fp) ? 0 : -1;
    fclose(fp);
    return k;
}

#endif
