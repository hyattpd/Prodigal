/*******************************************************************************
    fuzz_find_genes: libFuzzer target for the sequence path. The first byte
    of an input picks one of four contexts: metagenomic mode or single mode
    from ref_train.bin, each plain or with closed ends and masking. The rest
    is the sequence, given to prodigal_set_sequence and searched with
    prodigal_find_genes; the genes must lie in the sequence with every field
    in range, and their proteins must pass api_proteins_ok (api_gff.h). Run
    from the repository root (make fuzz; docs/testing.md).

    The contexts serve every input, as a GPL-boundary worker's does, so a
    failure can depend on the inputs before it and not reproduce alone.
*******************************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prodigal.h"
#include "api_gff.h"
#include "fasta.h"

#define MODEL_PATH "testdata/ground_truth/ref_train.bin"

static prodigal_ctx_t *ctx[4];

static void die(const char *what) {
    fprintf(stderr, "fuzz_find_genes: %s\n", what);
    exit(1);
}

int LLVMFuzzerInitialize(int *argc, char ***argv) {
    size_t len;
    void *model = load_file(MODEL_PATH, &len);
    int i;
    (void)argc; (void)argv;
    if (model == NULL)
        die("cannot read " MODEL_PATH " (run from the repository root)");
    for (i = 0; i < 4; i++) {
        prodigal_config_t cfg;
        prodigal_config_init(&cfg);
        cfg.meta_mode = i & 1;
        cfg.closed_ends = cfg.mask_regions = i >> 1;
        cfg.want_proteins = 1;
        if ((ctx[i] = prodigal_create(&cfg)) == NULL) die("prodigal_create");
        if (!cfg.meta_mode &&
            prodigal_load_training(ctx[i], model, len) != PRODIGAL_OK)
            die(prodigal_last_error(ctx[i]));
    }
    free(model);
    return 0;
}

/* Abort, so libFuzzer records the input, unless cond holds. */
#define CHECK(cond) do { if (!(cond)) { \
    fprintf(stderr, "fuzz_find_genes: %s fails\n", #cond); abort(); } } while (0)

/* Every gene lies in the sequence, every field is in its range, and the
   proteins are well formed. */
static void check_genes(prodigal_ctx_t *c, const prodigal_genes_soa_t *g,
                        const prodigal_stats_t *st) {
    prodigal_seq_info_t info;
    char why[128];
    int32_t i;
    prodigal_get_seq_info(c, &info);
    CHECK(g->n_genes >= 0 && st->n_genes == g->n_genes);
    CHECK(memchr(st->best_meta_desc, '\0', sizeof(st->best_meta_desc)) != NULL);
    for (i = 0; i < g->n_genes; i++) {
        CHECK(1 <= g->begin[i] && g->begin[i] <= g->end[i] &&
              g->end[i] <= info.length);
        CHECK((g->end[i] - g->begin[i] + 1) % 3 == 0);
        CHECK(g->strand[i] == 1 || g->strand[i] == -1);
        CHECK((uint32_t)g->partial_left[i] <= 1 && (uint32_t)g->partial_right[i] <= 1);
        CHECK(0 <= g->start_type[i] && g->start_type[i] <= 3);
        CHECK(50.0 <= g->confidence[i] && g->confidence[i] <= 100.0);
        CHECK(g->rbs_motif[i] != NULL && g->rbs_spacer[i] != NULL);
    }
    if (!api_proteins_ok(g, why, sizeof(why))) {
        fprintf(stderr, "fuzz_find_genes: %s\n", why);
        abort();
    }
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    prodigal_genes_soa_t *genes = NULL;
    prodigal_stats_t stats;
    prodigal_ctx_t *c;
    if (size < 1 || size - 1 > INT32_MAX) return 0;
    c = ctx[data[0] & 3];
    if (prodigal_set_sequence(c, (const char *)data + 1, (int32_t)(size - 1), NULL)
            == PRODIGAL_OK &&
        prodigal_find_genes(c, &genes, &stats) == PRODIGAL_OK) {
        check_genes(c, genes, &stats);
        prodigal_genes_free(genes);
    }
    return 0;
}
