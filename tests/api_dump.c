/*******************************************************************************
    api_dump: run the public library API over a FASTA file and print, per
    contig, a model line and one CLI-format GFF CDS line per gene
    (tests/api_gff.h), so the library can be diffed against the CLI
    (tests/api_compare.sh).

    usage: api_dump [--isolated] [-p single|meta] [-c] [-m] [-n] [-g TABLE]
                    [-t TRAINING] FASTA
           api_dump --from-cli [--isolated] [-p single|meta] GFF

    The options mean what they mean to the CLI. Single mode without -t trains
    on all records, as the CLI does. --from-cli instead translates a CLI GFF
    file into the same form. --isolated neutralises the known defects A and
    C (docs/testing.md): training sets the GC by hand, and every motif prints
    as "?". Only the public API (prodigal.h) is used.
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prodigal.h"
#include "api_gff.h"
#include "fasta.h"

static prodigal_ctx_t *ctx;

static void die(const char *what) {
    fprintf(stderr, "api_dump: %s%s%s\n", what, ctx ? ": " : "",
            ctx ? prodigal_last_error(ctx) : "");
    exit(1);
}

/* arg is non-NULL under --isolated. */
static void put_line(const char *line, void *arg) {
    char buf[API_GFF_LINE_MAX];
    if (arg != NULL) { api_mask_motif(strcpy(buf, line)); line = buf; }
    puts(line);
}

static void load_model(const char *path) {
    size_t n;
    void *buf = load_file(path, &n);
    if (buf == NULL) die("cannot read training file");
    if (prodigal_load_training(ctx, buf, n) != PRODIGAL_OK) die("load_training");
    free(buf);
}

int main(int argc, char **argv) {
    prodigal_config_t cfg;
    const char *train_path = NULL, *in = NULL;
    int meta = 0, cli = 0, isolated = 0, i;
    fasta_t fa;

    prodigal_config_init(&cfg);
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--from-cli") == 0) cli = 1;
        else if (strcmp(argv[i], "--isolated") == 0) isolated = 1;
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc)
            meta = strcmp(argv[++i], "meta") == 0;
        else if (strcmp(argv[i], "-c") == 0) cfg.closed_ends = 1;
        else if (strcmp(argv[i], "-m") == 0) cfg.mask_regions = 1;
        else if (strcmp(argv[i], "-n") == 0) cfg.force_nonsd = 1;
        else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc)
            cfg.trans_table = atoi(argv[++i]);
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc)
            train_path = argv[++i];
        else if (argv[i][0] != '-' && in == NULL) in = argv[i];
        else { fprintf(stderr, "api_dump: bad argument %s\n", argv[i]); return 2; }
    }
    if (in == NULL) { fprintf(stderr, "api_dump: no input\n"); return 2; }
    if (cli) {
        if (api_emit_cli_file(in, meta, put_line, isolated ? &isolated : NULL) != 0)
            die("cannot open GFF");
        return 0;
    }

    if (fasta_load(in, &fa) != 0) die("cannot read FASTA");
    cfg.meta_mode = meta;
    if ((ctx = prodigal_create(&cfg)) == NULL) die("prodigal_create failed");
    if (!meta) {
        if (train_path != NULL) load_model(train_path);
        else if (api_train(ctx, &fa, isolated) != PRODIGAL_OK) die("train");
    }

    if (api_emit_records(ctx, &fa, 0, fa.n, meta, put_line,
                         isolated ? &isolated : NULL) != PRODIGAL_OK)
        die("find_genes");
    fasta_free(&fa);
    prodigal_destroy(ctx);
    return 0;
}
