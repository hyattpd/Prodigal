/*******************************************************************************
    api_dump: run the public library API over a FASTA file and print, per
    contig, a model line and one CLI-format GFF CDS line per gene
    (tests/api_gff.h), so the library can be diffed against the CLI
    (tests/api_compare.sh).

    usage: api_dump [-p single|meta] [-c] [-m] [-n] [-g TABLE] [-t TRAINING]
                    [--write-model FILE] [--proteins] FASTA
           api_dump --from-cli [-p single|meta] [--faa FAA] GFF

    The options mean what they mean to the CLI. Single mode without -t trains
    on all records, as the CLI does; --write-model also writes that model to
    FILE, which then equals the CLI's -t file. --proteins adds each gene's
    protein (want_proteins). --from-cli instead translates a CLI GFF file
    into the same form, with the proteins of the CLI's -a file FAA. Only the
    public API (prodigal.h) is used.
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

static void put_line(const char *line, void *arg) { (void)arg; puts(line); }

static void load_model(const char *path) {
    size_t n;
    void *buf = load_file(path, &n);
    if (buf == NULL) die("cannot read training file");
    if (prodigal_load_training(ctx, buf, n) != PRODIGAL_OK) die("load_training");
    free(buf);
}

static void write_model(const char *path) {
    size_t n;
    void *buf;
    FILE *f;
    if (prodigal_export_training(ctx, &buf, &n) != PRODIGAL_OK) die("export_training");
    if ((f = fopen(path, "wb")) == NULL || fwrite(buf, 1, n, f) != n ||
        fclose(f) != 0)
        die("cannot write the model");
    free(buf);
}

int main(int argc, char **argv) {
    prodigal_config_t cfg;
    const char *train_path = NULL, *model_out = NULL, *faa = NULL, *in = NULL;
    int meta = 0, cli = 0, i;
    fasta_t fa;

    prodigal_config_init(&cfg);
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--from-cli") == 0) cli = 1;
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc)
            meta = strcmp(argv[++i], "meta") == 0;
        else if (strcmp(argv[i], "-c") == 0) cfg.closed_ends = 1;
        else if (strcmp(argv[i], "-m") == 0) cfg.mask_regions = 1;
        else if (strcmp(argv[i], "-n") == 0) cfg.force_nonsd = 1;
        else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc)
            cfg.trans_table = atoi(argv[++i]);
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc)
            train_path = argv[++i];
        else if (strcmp(argv[i], "--write-model") == 0 && i + 1 < argc)
            model_out = argv[++i];
        else if (strcmp(argv[i], "--proteins") == 0) cfg.want_proteins = 1;
        else if (strcmp(argv[i], "--faa") == 0 && i + 1 < argc) faa = argv[++i];
        else if (argv[i][0] != '-' && in == NULL) in = argv[i];
        else { fprintf(stderr, "api_dump: bad argument %s\n", argv[i]); return 2; }
    }
    if (in == NULL) { fprintf(stderr, "api_dump: no input\n"); return 2; }
    if (cli) {
        int rc = api_emit_cli_file(in, faa, meta, put_line, NULL);
        if (rc == API_CLI_MISMATCH) die("the GFF and -a files hold different genes");
        if (rc != 0) die("cannot read the GFF or -a file");
        return 0;
    }

    if (fasta_load(in, &fa) != 0) die("cannot read FASTA");
    cfg.meta_mode = meta;
    if ((ctx = prodigal_create(&cfg)) == NULL) die("prodigal_create failed");
    if (!meta) {
        if (train_path != NULL) load_model(train_path);
        else if (api_train(ctx, &fa) != PRODIGAL_OK) die("train");
        else if (model_out != NULL) write_model(model_out);
    }

    if (api_emit_records(ctx, &fa, 0, fa.n, meta, put_line, NULL) != PRODIGAL_OK)
        die("find_genes");
    fasta_free(&fa);
    prodigal_destroy(ctx);
    return 0;
}
