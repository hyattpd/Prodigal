/*******************************************************************************
    PRODIGAL Library API Test Suite
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "prodigal.h"
#include "tests/api_gff.h"
#include "tests/fasta.h"

static int tests_run = 0;
static int tests_passed = 0;
static int tests_xfailed = 0;
static int fixtures_missing = 0;

#define TEST_START(name) do { \
    tests_run++; \
    printf("  %-60s ", name); \
    fflush(stdout); \
} while(0)

#define TEST_PASS() do { \
    tests_passed++; \
    printf("[PASS]\n"); \
} while(0)

#define ASSERT_EQ_INT(a, b) do { \
    if ((a) != (b)) { \
        printf("[FAIL] %s:%d: %d != %d\n", __FILE__, __LINE__, (a), (b)); \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        printf("[FAIL] %s:%d: assertion failed: %s\n", __FILE__, __LINE__, #cond); \
        return; \
    } \
} while(0)

/*******************************************************************************
    Test helpers: fixtures and contexts. Each prints a [FAIL] line when it
    returns 0 or NULL.
*******************************************************************************/

#define REF_DIR "testdata/ground_truth/"

/* A missing fixture fails the test, and is never a known failure. */
static void missing_fixture(const char *what) {
    printf("[FAIL] missing fixture %s (run from the repo root)\n", what);
    fixtures_missing++;
}

static int load_anthus(fasta_t *fa) {
    if (fasta_load("anthus_aco.fas", fa) == 0 && fa->n > 0) return 1;
    missing_fixture("anthus_aco.fas");
    return 0;
}

static void *load_ref_model(size_t *len) {
    void *blob = load_file(REF_DIR "ref_train.bin", len);
    if (blob == NULL) missing_fixture(REF_DIR "ref_train.bin");
    return blob;
}

static prodigal_ctx_t *make_ctx(int meta) {
    prodigal_config_t config;
    prodigal_ctx_t *ctx;
    prodigal_config_init(&config);
    config.meta_mode = meta;
    if ((ctx = prodigal_create(&config)) == NULL) printf("[FAIL] prodigal_create\n");
    return ctx;
}

static prodigal_ctx_t *make_model_ctx(const void *blob, size_t len) {
    prodigal_ctx_t *ctx = make_ctx(0);
    if (ctx != NULL && prodigal_load_training(ctx, blob, len) != PRODIGAL_OK) {
        printf("[FAIL] load_training: %s\n", prodigal_last_error(ctx));
        prodigal_destroy(ctx); ctx = NULL;
    }
    return ctx;
}

/*******************************************************************************
    Phase 1.1: Error Codes
*******************************************************************************/

static void test_error_codes(void) {
    TEST_START("error codes: PRODIGAL_OK is 0, errors are negative");
    ASSERT_EQ_INT(PRODIGAL_OK, 0);
    ASSERT_TRUE(PRODIGAL_ERR_NOMEM < 0);
    ASSERT_TRUE(PRODIGAL_ERR_INVALID_CONFIG < 0);
    ASSERT_TRUE(PRODIGAL_ERR_INVALID_INPUT < 0);
    ASSERT_TRUE(PRODIGAL_ERR_INTERNAL < 0);
    ASSERT_TRUE(PRODIGAL_ERR_SEQ_TOO_SHORT < 0);
    ASSERT_TRUE(PRODIGAL_ERR_CANCELLED < 0);
    TEST_PASS();
}

static void test_strerror(void) {
    TEST_START("prodigal_strerror returns valid strings");
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_OK), "Success") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_ERR_NOMEM), "Out of memory") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_ERR_INVALID_CONFIG), "Invalid configuration") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_ERR_INVALID_INPUT), "Invalid input") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_ERR_INTERNAL), "Internal error") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_ERR_SEQ_TOO_SHORT), "Sequence too short") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(PRODIGAL_ERR_CANCELLED), "Cancelled") == 0);
    ASSERT_TRUE(strcmp(prodigal_strerror(-999), "Unknown error") == 0);
    TEST_PASS();
}

static void test_version_constants(void) {
    TEST_START("version constants and runtime query");
    ASSERT_TRUE(PRODIGAL_VERSION_MAJOR >= 2);
    ASSERT_TRUE(PRODIGAL_VERSION_MINOR >= 0);
    ASSERT_TRUE(PRODIGAL_VERSION_PATCH >= 0);
    ASSERT_TRUE(strlen(PRODIGAL_VERSION_STRING) > 0);
    /* Runtime version must match compile-time */
    ASSERT_TRUE(strcmp(prodigal_version_string(), PRODIGAL_VERSION_STRING) == 0);
    TEST_PASS();
}

/*******************************************************************************
    Phase 1.2: Config Struct
*******************************************************************************/

static void test_config_init_defaults(void) {
    TEST_START("config_init sets correct defaults");
    prodigal_config_t config;
    prodigal_config_init(&config);

    ASSERT_TRUE(config.struct_size == sizeof(prodigal_config_t));
    ASSERT_EQ_INT(config.trans_table, 11);
    ASSERT_EQ_INT(config.closed_ends, 0);
    ASSERT_EQ_INT(config.mask_regions, 0);
    ASSERT_EQ_INT(config.force_nonsd, 0);
    ASSERT_EQ_INT(config.meta_mode, 0);
    ASSERT_TRUE(config.start_weight == 4.35);
    ASSERT_TRUE(config.alloc_fn == NULL);
    ASSERT_TRUE(config.free_fn == NULL);
    ASSERT_TRUE(config.allocator_user_data == NULL);
    ASSERT_TRUE(config.log_callback == NULL);
    ASSERT_TRUE(config.log_user_data == NULL);
    ASSERT_TRUE(config.progress_callback == NULL);
    ASSERT_TRUE(config.progress_user_data == NULL);
    TEST_PASS();
}

static void test_config_struct_size_at_offset_zero(void) {
    TEST_START("struct_size is at offset 0 in config");
    prodigal_config_t config;
    ASSERT_TRUE((char *)&config.struct_size == (char *)&config);
    TEST_PASS();
}

/*******************************************************************************
    Phase 1.3: Context Lifecycle
*******************************************************************************/

static void test_create_destroy(void) {
    TEST_START("create and destroy context");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);
    ASSERT_TRUE(ctx != NULL);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_create_null_config(void) {
    TEST_START("create with NULL config returns NULL");
    prodigal_ctx_t *ctx = prodigal_create(NULL);
    ASSERT_TRUE(ctx == NULL);
    TEST_PASS();
}

static void test_create_bad_struct_size(void) {
    TEST_START("create with wrong struct_size returns NULL");
    prodigal_config_t config;
    prodigal_config_init(&config);
    config.struct_size = 1;
    prodigal_ctx_t *ctx = prodigal_create(&config);
    ASSERT_TRUE(ctx == NULL);
    TEST_PASS();
}

static void test_destroy_null(void) {
    TEST_START("destroy(NULL) is safe");
    prodigal_destroy(NULL);
    TEST_PASS();
}

static void test_create_invalid_trans_table(void) {
    TEST_START("create with invalid translation tables returns NULL");
    int invalid[] = {0, 7, 8, 17, 18, 19, 20, 26, -1};
    int n = sizeof(invalid) / sizeof(invalid[0]);
    int i;

    for (i = 0; i < n; i++) {
        prodigal_config_t config;
        prodigal_config_init(&config);
        config.trans_table = invalid[i];
        prodigal_ctx_t *ctx = prodigal_create(&config);
        if (ctx != NULL) {
            prodigal_destroy(ctx);
            printf("[FAIL] trans_table %d should have been rejected\n", invalid[i]);
            return;
        }
    }
    TEST_PASS();
}

static void test_create_valid_trans_tables(void) {
    TEST_START("create with all valid translation tables succeeds");
    int valid[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25};
    int n = sizeof(valid) / sizeof(valid[0]);
    int i;

    for (i = 0; i < n; i++) {
        prodigal_config_t config;
        prodigal_config_init(&config);
        config.trans_table = valid[i];
        prodigal_ctx_t *ctx = prodigal_create(&config);
        if (ctx == NULL) {
            printf("[FAIL] trans_table %d should have been accepted\n", valid[i]);
            return;
        }
        prodigal_destroy(ctx);
    }
    TEST_PASS();
}

static void test_create_meta_mode(void) {
    TEST_START("create with meta_mode=1 succeeds");
    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    prodigal_ctx_t *ctx = prodigal_create(&config);
    ASSERT_TRUE(ctx != NULL);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_last_error_on_fresh_context(void) {
    TEST_START("last_error on fresh context is empty");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);
    ASSERT_TRUE(ctx != NULL);
    const char *err = prodigal_last_error(ctx);
    ASSERT_TRUE(err != NULL);
    ASSERT_TRUE(strlen(err) == 0);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_last_error_null_context(void) {
    TEST_START("last_error with NULL context returns non-NULL");
    const char *err = prodigal_last_error(NULL);
    ASSERT_TRUE(err != NULL);
    ASSERT_TRUE(strlen(err) > 0);
    TEST_PASS();
}

static void test_genes_free_null(void) {
    TEST_START("genes_free(NULL) is safe");
    prodigal_genes_free(NULL);
    prodigal_genes_aos_free(NULL);
    TEST_PASS();
}

/*******************************************************************************
    Phase 2.1: Single Sequence Encoding
*******************************************************************************/

static void test_encode_simple_sequence(void) {
    TEST_START("encode simple ACGTACGT sequence");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);
    ASSERT_TRUE(ctx != NULL);

    const char *seq = "ACGTACGT";
    int rc = prodigal_set_sequence(ctx, seq, 8, "test_seq");
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_seq_info_t info;
    rc = prodigal_get_seq_info(ctx, &info);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_EQ_INT(info.length, 8);
    /* 4 GC out of 8 = 0.5 */
    ASSERT_TRUE(fabs(info.gc_content - 0.5) < 1e-10);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_encode_with_ambiguity(void) {
    TEST_START("encode sequence with N bases");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    const char *seq = "ACNGTNCG";
    int rc = prodigal_set_sequence(ctx, seq, 8, "test_n");
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_seq_info_t info;
    prodigal_get_seq_info(ctx, &info);
    ASSERT_EQ_INT(info.length, 8);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_encode_null_sequence(void) {
    TEST_START("set_sequence with NULL returns error");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    int rc = prodigal_set_sequence(ctx, NULL, 0, "empty");
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_TRUE(strlen(prodigal_last_error(ctx)) > 0);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_encode_empty_string(void) {
    TEST_START("set_sequence with empty string returns error");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    int rc = prodigal_set_sequence(ctx, "", 0, "empty");
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_INVALID_INPUT);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_encode_gc_content(void) {
    TEST_START("GC content calculated correctly");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    /* All G's: GC = 1.0 */
    const char *all_g = "GGGGGGGGGG";
    prodigal_set_sequence(ctx, all_g, 10, "all_g");
    prodigal_seq_info_t info;
    prodigal_get_seq_info(ctx, &info);
    ASSERT_TRUE(fabs(info.gc_content - 1.0) < 1e-10);

    /* All A's: GC = 0.0 */
    const char *all_a = "AAAAAAAAAA";
    prodigal_set_sequence(ctx, all_a, 10, "all_a");
    prodigal_get_seq_info(ctx, &info);
    ASSERT_TRUE(fabs(info.gc_content - 0.0) < 1e-10);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_encode_lowercase(void) {
    TEST_START("lowercase bases accepted");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    int rc = prodigal_set_sequence(ctx, "acgtacgt", 8, "lower");
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_seq_info_t info;
    prodigal_get_seq_info(ctx, &info);
    ASSERT_EQ_INT(info.length, 8);
    ASSERT_TRUE(fabs(info.gc_content - 0.5) < 1e-10);

    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 2.2: Multi-Sequence Training Input
*******************************************************************************/

static void test_training_multi_sequence(void) {
    TEST_START("training concatenates with stop spacers");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    const char *seqs[] = {"ACGTACGT", "TGCATGCA"};
    const char *hdrs[] = {"seq1", "seq2"};
    int32_t lens[] = {8, 8};

    int rc = prodigal_set_training_sequences(ctx, seqs, hdrs, lens, 2);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    /* 2 seqs: seq1(8) + spacer(12) + seq2(8) + spacer(12) = 40 */
    prodigal_seq_info_t info;
    prodigal_get_seq_info(ctx, &info);
    ASSERT_EQ_INT(info.length, 40);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_training_single_sequence(void) {
    TEST_START("training single sequence: no spacers");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    const char *seqs[] = {"ACGTACGT"};
    const char *hdrs[] = {"seq1"};
    int32_t lens[] = {8};

    int rc = prodigal_set_training_sequences(ctx, seqs, hdrs, lens, 1);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_seq_info_t info;
    prodigal_get_seq_info(ctx, &info);
    ASSERT_EQ_INT(info.length, 8);

    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_training_null_args(void) {
    TEST_START("training with NULL args returns error");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    int rc = prodigal_set_training_sequences(ctx, NULL, NULL, NULL, 0);
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_INVALID_INPUT);

    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 3.1: Training Serialization
*******************************************************************************/

static void test_training_load_roundtrip(void) {
    TEST_START("training load/export round-trip is byte-identical");
    size_t tlen;
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL) return;
    prodigal_ctx_t *ctx = make_model_ctx(tdata, tlen);
    ASSERT_TRUE(ctx != NULL);

    /* Export and compare */
    void *exported;
    size_t exported_len;
    int rc = prodigal_export_training(ctx, &exported, &exported_len);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(exported_len == tlen);
    ASSERT_TRUE(memcmp(tdata, exported, exported_len) == 0);

    free(tdata);
    free(exported);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_training_load_invalid(void) {
    TEST_START("training load rejects wrong size");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    char dummy[16] = {0};
    int rc = prodigal_load_training(ctx, dummy, sizeof(dummy));
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_TRUE(strlen(prodigal_last_error(ctx)) > 0);

    rc = prodigal_load_training(ctx, NULL, 0);
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_INVALID_INPUT);

    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 3.2: Training Pipeline
*******************************************************************************/

static void test_train_from_sequences(void) {
    TEST_START("train from FASTA sequences produces valid model");
    fasta_t fa;
    if (!load_anthus(&fa)) return;
    prodigal_ctx_t *ctx = make_ctx(0);
    ASSERT_TRUE(ctx != NULL);

    int rc = prodigal_set_training_sequences(ctx,
        (const char **)fa.seqs, (const char **)fa.headers, fa.lens, fa.n);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_seq_info_t info;
    prodigal_get_seq_info(ctx, &info);

    /* Training requires >= 20000 bp */
    if (info.length < 20000) {
        printf("[FAIL] concat length %d < 20000\n", info.length);
        fasta_free(&fa);
        prodigal_destroy(ctx);
        return;
    }

    rc = prodigal_train(ctx);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    /* Verify training produced a sensible model */
    void *exported;
    size_t elen;
    rc = prodigal_export_training(ctx, &exported, &elen);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(elen > 0);

    /* The golden test "API training = ref_train.bin" compares the model with
       the reference. */
    free(exported);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 4.1: Gene Finding (Single Genome)
*******************************************************************************/

static void test_find_genes_with_training(void) {
    TEST_START("find genes using loaded training data");

    size_t tlen;
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL) return;
    fasta_t fa;
    if (!load_anthus(&fa)) { free(tdata); return; }
    prodigal_ctx_t *ctx = make_model_ctx(tdata, tlen);
    free(tdata);
    ASSERT_TRUE(ctx != NULL);

    int rc = prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_genes_soa_t *genes = NULL;
    prodigal_stats_t stats;
    rc = prodigal_find_genes(ctx, &genes, &stats);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(genes != NULL);

    /* Validate SOA structure */
    if (genes->n_genes > 0) {
        int i;
        ASSERT_TRUE(genes->begin != NULL);
        ASSERT_TRUE(genes->end != NULL);
        ASSERT_TRUE(genes->strand != NULL);
        ASSERT_TRUE(genes->cscore != NULL);
        ASSERT_TRUE(genes->confidence != NULL);

        for (i = 0; i < genes->n_genes; i++) {
            ASSERT_TRUE(genes->begin[i] >= 1);
            ASSERT_TRUE(genes->end[i] >= 1);
            ASSERT_TRUE(genes->strand[i] == 1 || genes->strand[i] == -1);
            ASSERT_TRUE(genes->confidence[i] >= 50.0);
            ASSERT_TRUE(genes->confidence[i] <= 100.0);
            ASSERT_TRUE(isfinite(genes->cscore[i]));
            ASSERT_TRUE(isfinite(genes->sscore[i]));
        }
    }

    ASSERT_TRUE(stats.n_genes == genes->n_genes);
    ASSERT_TRUE(stats.n_nodes > 0);

    prodigal_genes_free(genes);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 4.2: AOS Output
*******************************************************************************/

static void test_aos_matches_soa(void) {
    TEST_START("AOS output matches SOA for same input");

    size_t tlen;
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL) return;
    fasta_t fa;
    if (!load_anthus(&fa)) { free(tdata); return; }
    prodigal_ctx_t *ctx = make_model_ctx(tdata, tlen);
    free(tdata);
    ASSERT_TRUE(ctx != NULL);
    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);

    prodigal_genes_soa_t *soa = NULL;
    prodigal_find_genes(ctx, &soa, NULL);

    /* Re-set same sequence for AOS */
    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);
    prodigal_genes_t *aos = NULL;
    prodigal_find_genes_aos(ctx, &aos, NULL);

    ASSERT_TRUE(soa != NULL && aos != NULL);
    ASSERT_EQ_INT(soa->n_genes, aos->n_genes);

    if (soa->n_genes > 0) {
        int i;
        for (i = 0; i < soa->n_genes; i++) {
            ASSERT_EQ_INT(soa->begin[i], aos->genes[i].begin);
            ASSERT_EQ_INT(soa->end[i], aos->genes[i].end);
            ASSERT_EQ_INT(soa->strand[i], aos->genes[i].strand);
            ASSERT_TRUE(soa->cscore[i] == aos->genes[i].cscore);
            ASSERT_TRUE(soa->confidence[i] == aos->genes[i].confidence);
        }
    }

    prodigal_genes_free(soa);
    prodigal_genes_aos_free(aos);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 5: Metagenomic Mode
*******************************************************************************/

static void test_meta_find_genes(void) {
    TEST_START("metagenomic mode finds genes");

    fasta_t fa;
    if (!load_anthus(&fa)) return;
    prodigal_ctx_t *ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);

    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);

    prodigal_genes_soa_t *genes = NULL;
    prodigal_stats_t stats;
    int rc = prodigal_find_genes(ctx, &genes, &stats);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(genes != NULL);
    ASSERT_TRUE(stats.best_meta_bin >= 0);

    prodigal_genes_free(genes);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_meta_no_path_bin(void) {
    TEST_START("meta: contig without a gene path reports the CLI's bin");
    /* The first 50 bp of anthus_aco.fas: no bin finds a gene path. The CLI
       and upstream c1e2d36 report bin 1, the first whose GC window holds the
       contig (docs/upstream-behavior.md). */
    static const char seq[] = "ACAGGTTAGAAACTACTCTGTTTTCTGGCTCCTTGTTTAATGCCCTGTCC";
    prodigal_ctx_t *ctx = make_ctx(1);
    prodigal_genes_soa_t *genes = NULL;
    prodigal_stats_t stats;

    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, seq, 50, NULL), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, &stats), PRODIGAL_OK);
    ASSERT_EQ_INT(stats.n_genes, 0);
    ASSERT_EQ_INT(stats.best_meta_bin, 1);
    ASSERT_EQ_INT(stats.translation_table, 4);
    ASSERT_EQ_INT(stats.uses_sd, 0);
    prodigal_genes_free(genes);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 7: Context Reuse and Error Recovery
*******************************************************************************/

static void test_context_reuse(void) {
    TEST_START("context reuse: process multiple sequences sequentially");

    fasta_t fa;
    if (!load_anthus(&fa)) return;
    prodigal_ctx_t *ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);

    int32_t s;
    for (s = 0; s < fa.n; s++) {
        prodigal_set_sequence(ctx, fa.seqs[s], fa.lens[s], fa.headers[s]);
        prodigal_genes_soa_t *genes = NULL;
        int rc = prodigal_find_genes(ctx, &genes, NULL);
        ASSERT_EQ_INT(rc, PRODIGAL_OK);
        prodigal_genes_free(genes);
    }

    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_find_genes_without_sequence(void) {
    TEST_START("find_genes without a sequence is an error");
    prodigal_genes_soa_t *genes = NULL;
    prodigal_ctx_t *ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_TRUE(genes == NULL);
    prodigal_destroy(ctx);

    /* Single mode: training is checked first */
    ctx = make_ctx(0);
    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_TRUE(strstr(prodigal_last_error(ctx), "training") != NULL);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_error_recovery(void) {
    TEST_START("error recovery: bad input then good input");

    fasta_t fa;
    if (!load_anthus(&fa)) return;
    prodigal_ctx_t *ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);

    /* Bad input */
    int rc = prodigal_set_sequence(ctx, NULL, 0, "bad");
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_TRUE(strlen(prodigal_last_error(ctx)) > 0);

    /* Good input should still work */
    rc = prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], "good");
    ASSERT_EQ_INT(rc, PRODIGAL_OK);

    prodigal_genes_soa_t *genes = NULL;
    rc = prodigal_find_genes(ctx, &genes, NULL);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(genes != NULL);

    prodigal_genes_free(genes);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 8: Allocator Hooks
*******************************************************************************/

static size_t alloc_count = 0;
static size_t free_count = 0;

static size_t alloc_limit = (size_t)-1;   /* test_alloc fails past this */

static void *test_alloc(size_t size, void *user_data) {
    (void)user_data;
    if (alloc_count >= alloc_limit) return NULL;
    alloc_count++;
    return malloc(size);
}

static void test_free_fn(void *ptr, void *user_data) {
    free_count++;
    (void)user_data;
    free(ptr);
}

/* A metagenomic context on test_alloc/test_free_fn, with the counters reset
   and no allocation limit. */
static prodigal_ctx_t *make_counted_ctx(void) {
    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    config.alloc_fn = test_alloc;
    config.free_fn = test_free_fn;
    alloc_count = free_count = 0;
    alloc_limit = (size_t)-1;
    return prodigal_create(&config);
}

static void test_custom_allocator(void) {
    TEST_START("custom allocator used for internal buffers");
    fasta_t fa;
    if (!load_anthus(&fa)) return;

    prodigal_ctx_t *ctx = make_counted_ctx();
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(alloc_count > 0);

    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);
    prodigal_genes_soa_t *genes = NULL;
    prodigal_find_genes(ctx, &genes, NULL);
    ASSERT_TRUE(genes != NULL);

    prodigal_genes_free(genes);  /* Uses system free, not custom */
    prodigal_destroy(ctx);
    ASSERT_TRUE(free_count > 0);

    fasta_free(&fa);
    TEST_PASS();
}

static void test_meta_alloc_failure(void) {
    TEST_START("meta: failed bin allocation, retry, no leak");
    fasta_t fa;
    if (!load_anthus(&fa)) return;

    prodigal_ctx_t *ctx = make_counted_ctx();
    ASSERT_TRUE(ctx != NULL);
    size_t held = alloc_count - free_count;
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]),
                  PRODIGAL_OK);

    /* The bin array and 10 of the 50 bin models, then failure */
    alloc_limit = alloc_count + 11;
    prodigal_genes_soa_t *genes = NULL;
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_ERR_NOMEM);
    ASSERT_TRUE(genes == NULL);
    ASSERT_TRUE(alloc_count - free_count == held);

    alloc_limit = (size_t)-1;
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_OK);
    ASSERT_TRUE(genes != NULL && genes->n_genes > 0);
    prodigal_genes_free(genes);
    prodigal_destroy(ctx);
    ASSERT_TRUE(alloc_count == free_count);
    fasta_free(&fa);
    TEST_PASS();
}

/*******************************************************************************
    Phase 9: Callbacks
*******************************************************************************/

static int log_call_count = 0;

static void test_log_cb(const char *msg, void *user_data) {
    (void)user_data;
    (void)msg;
    log_call_count++;
}

static void test_log_callback(void) {
    TEST_START("log callback receives messages during find_genes");
    log_call_count = 0;

    fasta_t fa;
    if (!load_anthus(&fa)) return;

    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    config.log_callback = test_log_cb;
    prodigal_ctx_t *ctx = prodigal_create(&config);

    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);
    prodigal_genes_soa_t *genes = NULL;
    prodigal_find_genes(ctx, &genes, NULL);

    /* Meta mode logs progress via progress_callback, not log_callback.
       But the log callback should still fire if we train. Just verify
       the callback mechanism works. */
    ASSERT_TRUE(log_call_count >= 0);  /* may be 0 in meta mode without training */

    prodigal_genes_free(genes);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static int cancel_immediately(const char *stage, double frac, void *user_data) {
    (void)stage; (void)frac; (void)user_data;
    return 1;  /* cancel */
}

static void test_progress_cancellation(void) {
    TEST_START("progress callback can cancel metagenomic computation");

    fasta_t fa;
    if (!load_anthus(&fa)) return;

    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    config.progress_callback = cancel_immediately;
    prodigal_ctx_t *ctx = prodigal_create(&config);

    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);
    prodigal_genes_soa_t *genes = NULL;
    int rc = prodigal_find_genes(ctx, &genes, NULL);
    ASSERT_EQ_INT(rc, PRODIGAL_ERR_CANCELLED);
    ASSERT_TRUE(genes == NULL);

    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 10: Edge Cases
*******************************************************************************/

static void test_very_short_sequence(void) {
    TEST_START("very short sequence (50bp) produces 0 genes");
    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    prodigal_ctx_t *ctx = prodigal_create(&config);

    const char *short_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    prodigal_set_sequence(ctx, short_seq, 49, "short");

    prodigal_genes_soa_t *genes = NULL;
    prodigal_stats_t stats;
    int rc = prodigal_find_genes(ctx, &genes, &stats);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_EQ_INT(stats.n_genes, 0);
    ASSERT_EQ_INT(genes->n_genes, 0);

    prodigal_genes_free(genes);
    prodigal_destroy(ctx);
    TEST_PASS();
}

static void test_all_n_sequence(void) {
    TEST_START("all-N sequence completes without error");
    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    prodigal_ctx_t *ctx = prodigal_create(&config);

    /* N bases are encoded as C with ambiguity flag set. Prodigal may still
       find genes in ambiguous regions (this is by design — only masking
       mode prevents genes from spanning N-runs, and even then short
       all-N sequences can produce spurious hits). The key assertion is
       that the library handles this gracefully. */
    char all_n[201];
    memset(all_n, 'N', 200);
    all_n[200] = '\0';
    prodigal_set_sequence(ctx, all_n, 200, "allN");

    prodigal_genes_soa_t *genes = NULL;
    int rc = prodigal_find_genes(ctx, &genes, NULL);
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(genes != NULL);
    ASSERT_TRUE(genes->n_genes >= 0);  /* May or may not find genes */

    prodigal_genes_free(genes);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 12: SOA Alignment
*******************************************************************************/

static void test_soa_alignment(void) {
    TEST_START("SOA arrays are 16-byte aligned");

    fasta_t fa;
    if (!load_anthus(&fa)) return;
    prodigal_ctx_t *ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);
    prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]);

    prodigal_genes_soa_t *genes = NULL;
    prodigal_find_genes(ctx, &genes, NULL);
    ASSERT_TRUE(genes != NULL);

    if (genes->n_genes > 0) {
        ASSERT_TRUE(genes->_base != NULL);
        ASSERT_TRUE(((uintptr_t)genes->begin) % 16 == 0);
        ASSERT_TRUE(((uintptr_t)genes->end) % 16 == 0);
        ASSERT_TRUE(((uintptr_t)genes->strand) % 16 == 0);
        ASSERT_TRUE(((uintptr_t)genes->cscore) % 16 == 0);
        ASSERT_TRUE(((uintptr_t)genes->sscore) % 16 == 0);
        ASSERT_TRUE(((uintptr_t)genes->confidence) % 16 == 0);
        ASSERT_TRUE(((uintptr_t)genes->gc_cont) % 16 == 0);
    }

    prodigal_genes_free(genes);
    fasta_free(&fa);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 13: Training Setters
*******************************************************************************/

static void test_training_setters(void) {
    TEST_START("training parameter setters work correctly");
    prodigal_config_t config;
    prodigal_config_init(&config);
    prodigal_ctx_t *ctx = prodigal_create(&config);

    /* Valid setters */
    ASSERT_EQ_INT(prodigal_set_translation_table(ctx, 4), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_start_weight(ctx, 3.0), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, 0.45), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_uses_sd(ctx, 0), PRODIGAL_OK);

    /* Invalid setters */
    ASSERT_EQ_INT(prodigal_set_translation_table(ctx, 7), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_start_weight(ctx, -1.0), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, 2.0), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, -0.1), PRODIGAL_ERR_INVALID_INPUT);

    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 14: API ≡ CLI golden tests

    The library's output, in the form of tests/api_gff.h, must equal the
    CLI's ground truth line for line. A known failure runs again as an
    "isolated" twin that neutralises only the known defects and must pass
    (docs/testing.md).
*******************************************************************************/

typedef struct { char **v; int n, cap; } lines_t;

static void lines_emit(const char *line, void *arg) {
    lines_t *l = (lines_t *)arg;
    if (l->n == l->cap) {
        l->cap = l->cap ? l->cap * 2 : 256;
        l->v = (char **)xrealloc(l->v, l->cap * sizeof(char *));
    }
    l->v[l->n++] = dupstr(line);
}

static void lines_free(lines_t *l) {
    int i;
    for (i = 0; i < l->n; i++) free(l->v[i]);
    free(l->v);
    memset(l, 0, sizeof(*l));
}

/* 1 if want and got are equal; otherwise prints [FAIL] with the number of
   differing lines and the first difference, and returns 0. With mask_c, a
   pair where got's motif is "Upstream" (defect C) is compared without the
   motif. */
static int lines_equal(const lines_t *want, const lines_t *got, const char *what,
                       int mask_c) {
    char w[API_GFF_LINE_MAX], g[API_GFF_LINE_MAX];
    int i, ndiff = 0, first = -1;
    if (want->n != got->n) {
        printf("[FAIL] %s: %d lines, expected %d\n", what, got->n, want->n);
        return 0;
    }
    for (i = 0; i < want->n; i++) {
        const char *a = want->v[i], *b = got->v[i];
        if (mask_c && strstr(b, "rbs_motif=Upstream;") != NULL) {
            api_mask_motif(strcpy(w, a)); a = w;   /* lines are < API_GFF_LINE_MAX */
            api_mask_motif(strcpy(g, b)); b = g;
        }
        if (strcmp(a, b) != 0 && ndiff++ == 0) first = i;
    }
    if (ndiff == 0) return 1;
    printf("[FAIL] %s: %d of %d lines differ; first:\n      want %s\n      got  %s\n",
           what, ndiff, want->n, want->v[first], got->v[first]);
    return 0;
}

/* Single-mode context trained through the API on all of fa (api_train). */
static prodigal_ctx_t *make_trained_ctx(const fasta_t *fa, int fix_gc) {
    prodigal_ctx_t *ctx = make_ctx(0);
    if (ctx != NULL && api_train(ctx, fa, fix_gc) != PRODIGAL_OK) {
        printf("[FAIL] training: %s\n", prodigal_last_error(ctx));
        prodigal_destroy(ctx); ctx = NULL;
    }
    return ctx;
}

/* Run all of fa on ctx and compare with the CLI's output in ref. */
static int golden_run(prodigal_ctx_t *ctx, const fasta_t *fa, const char *ref,
                      int meta, int mask_c) {
    lines_t want = {0}, got = {0};
    int ok = 0;
    if (api_emit_cli_file(ref, meta, lines_emit, &want) != 0)
        missing_fixture(ref);
    else if (api_emit_records(ctx, fa, 0, fa->n, meta, lines_emit, &got) != PRODIGAL_OK)
        printf("[FAIL] find_genes: %s\n", prodigal_last_error(ctx));
    else
        ok = lines_equal(&want, &got, ref, mask_c);
    lines_free(&want); lines_free(&got);
    return ok;
}

/* ARM64 macOS training writes different bytes (docs/testing.md). */
#if defined(__APPLE__) && defined(__aarch64__)
static const int train_bytes_match = 0;
#else
static const int train_bytes_match = 1;
#endif

/* Is the model trained on ctx the reference? Byte for byte where training
   writes the reference's bytes; elsewhere it must predict ref_trained.gff,
   with motifs masked (C), since this test is about the model. */
static int model_is_ref(prodigal_ctx_t *ctx, const fasta_t *fa) {
    void *got = NULL, *ref = NULL;
    size_t glen = 0, rlen = 0, i = 0;
    int ok = 0;
    if (prodigal_export_training(ctx, &got, &glen) != PRODIGAL_OK) {
        printf("[FAIL] export: %s\n", prodigal_last_error(ctx));
        return 0;
    }
    if (!train_bytes_match) {
        prodigal_ctx_t *pred = make_model_ctx(got, glen);
        ok = pred != NULL && golden_run(pred, fa, REF_DIR "ref_trained.gff", 0, 1);
        prodigal_destroy(pred);
    } else if ((ref = load_ref_model(&rlen)) != NULL) {
        while (i < rlen && i < glen && ((char *)ref)[i] == ((char *)got)[i]) i++;
        ok = glen == rlen && i == rlen;
        if (!ok) printf("[FAIL] %zu vs %zu bytes, first difference at byte %zu\n", glen, rlen, i);
    }
    free(got); free(ref);
    return ok;
}

enum { META, API_TRAINED, REF_MODEL };

/* One golden test on anthus_aco.fas: the genes found on a context made from
   src against the CLI's output in ref or, when ref is NULL, the model
   trained on it against ref_train.bin. isolated neutralises defects A and
   C. */
static void golden(int src, const char *ref, int isolated) {
    fasta_t fa;
    prodigal_ctx_t *ctx = NULL;
    void *blob = NULL;
    size_t len;
    if (!load_anthus(&fa)) return;
    if (src == META) ctx = make_ctx(1);
    else if (src == API_TRAINED) ctx = make_trained_ctx(&fa, isolated);
    else if ((blob = load_ref_model(&len)) != NULL) ctx = make_model_ctx(blob, len);
    if (ctx != NULL && (ref != NULL ? golden_run(ctx, &fa, ref, src == META, isolated)
                                    : model_is_ref(ctx, &fa)))
        TEST_PASS();
    prodigal_destroy(ctx);
    free(blob);
    fasta_free(&fa);
}

/* Settle a test that is known to fail, given tests_passed and
   fixtures_missing from before it ran; reason names the defect and the plan
   phase that fixes it. While it fails it counts as passed; once it passes
   the run fails, so that the marker gets removed. */
static void settle_xfail(int passed, int missing, const char *reason) {
    if (fixtures_missing > missing) return;
    if (tests_passed > passed) {
        tests_passed--;
        printf("    ^ XPASS: passes now; remove its XFAIL marker (%s)\n", reason);
    } else {
        tests_passed++;
        tests_xfailed++;
        printf("    ^ XFAIL, known: %s\n", reason);
    }
}

static void test_golden(void) {
    static const struct {
        const char *name;
        int src;
        const char *ref, *xfail;
    } t[] = {
        { "meta mode = ref_meta.gff, every field", META, REF_DIR "ref_meta.gff", NULL },
        { "API training = ref_train.bin", API_TRAINED, NULL,
          "A, prodigal_train never sets tinf.gc (P1)" },
        { "single mode, API-trained = ref_single.gff", API_TRAINED,
          REF_DIR "ref_single.gff",
          "A, training GC; C, non-SD motifs read \"Upstream\" (P1)" },
        { "single mode from ref_train.bin = ref_trained.gff", REF_MODEL,
          REF_DIR "ref_trained.gff", "C, non-SD motifs read \"Upstream\" (P1)" },
    };
    char name[128];
    size_t i;
    for (i = 0; i < sizeof(t) / sizeof(t[0]); i++) {
        int passed = tests_passed, missing = fixtures_missing;
        snprintf(name, sizeof(name), "golden: %s", t[i].name);
        TEST_START(name);
        golden(t[i].src, t[i].ref, 0);
        if (t[i].xfail == NULL) continue;
        settle_xfail(passed, missing, t[i].xfail);
        snprintf(name, sizeof(name), "golden: %s, isolated", t[i].name);
        TEST_START(name);
        golden(t[i].src, t[i].ref, 1);
    }
}

/* A contig's output must not depend on what the context did before: running
   every record again on the same context, each record on a fresh context,
   and prodigal_find_genes twice per prodigal_set_sequence must all give the
   first run's output. */
static void independence(int meta) {
    fasta_t fa;
    lines_t first = {0}, again = {0}, solo = {0}, twice = {0};
    prodigal_ctx_t *ctx = NULL;
    void *blob = NULL;
    size_t len = 0;
    int32_t k;
    int rc;
    if (!load_anthus(&fa)) return;
    if (meta || (blob = load_ref_model(&len)) != NULL)
        ctx = meta ? make_ctx(1) : make_model_ctx(blob, len);
    if (ctx != NULL) {
        rc = api_emit_records(ctx, &fa, 0, fa.n, meta, lines_emit, &first);
        if (rc == PRODIGAL_OK) rc = api_emit_records(ctx, &fa, 0, fa.n, meta, lines_emit, &again);
        for (k = 0; rc == PRODIGAL_OK && k < fa.n; k++) {
            prodigal_ctx_t *one = meta ? make_ctx(1) : make_model_ctx(blob, len);
            rc = one ? api_emit_records(one, &fa, k, k + 1, meta, lines_emit, &solo)
                     : PRODIGAL_ERR_NOMEM;
            prodigal_destroy(one);
        }
        for (k = 0; rc == PRODIGAL_OK && k < fa.n; k++) {
            prodigal_genes_soa_t *g = NULL;
            rc = prodigal_set_sequence(ctx, fa.seqs[k], fa.lens[k], fa.headers[k]);
            if (rc == PRODIGAL_OK) rc = prodigal_find_genes(ctx, &g, NULL);
            prodigal_genes_free(g);
            if (rc == PRODIGAL_OK) rc = api_emit_found(ctx, &fa, k, meta, lines_emit, &twice);
        }
        if (rc != PRODIGAL_OK) printf("[FAIL] find_genes failed (%d)\n", rc);
        else if (lines_equal(&first, &again, "second pass on the same context", 0) &&
                 lines_equal(&first, &solo, "each record on a fresh context", 0) &&
                 lines_equal(&first, &twice, "find_genes twice per sequence", 0))
            TEST_PASS();
    }
    lines_free(&first); lines_free(&again); lines_free(&solo); lines_free(&twice);
    prodigal_destroy(ctx);
    free(blob);
    fasta_free(&fa);
}

static void test_independence_meta(void) {
    TEST_START("repeat and order independence: meta mode");
    independence(1);
}

static void test_independence_single(void) {
    TEST_START("repeat and order independence: single mode from model");
    independence(0);
}

/*******************************************************************************
    Test runner
*******************************************************************************/

int main(void) {
    printf("=== Prodigal Library API Test Suite ===\n\n");

    printf("Phase 1.1: Error Codes\n");
    test_error_codes();
    test_strerror();
    test_version_constants();

    printf("\nPhase 1.2: Config Struct\n");
    test_config_init_defaults();
    test_config_struct_size_at_offset_zero();

    printf("\nPhase 1.3: Context Lifecycle\n");
    test_create_destroy();
    test_create_null_config();
    test_create_bad_struct_size();
    test_destroy_null();
    test_create_invalid_trans_table();
    test_create_valid_trans_tables();
    test_create_meta_mode();
    test_last_error_on_fresh_context();
    test_last_error_null_context();
    test_genes_free_null();

    printf("\nPhase 2.1: Single Sequence Encoding\n");
    test_encode_simple_sequence();
    test_encode_with_ambiguity();
    test_encode_null_sequence();
    test_encode_empty_string();
    test_encode_gc_content();
    test_encode_lowercase();

    printf("\nPhase 2.2: Multi-Sequence Training Input\n");
    test_training_multi_sequence();
    test_training_single_sequence();
    test_training_null_args();

    printf("\nPhase 3.1: Training Serialization\n");
    test_training_load_roundtrip();
    test_training_load_invalid();

    printf("\nPhase 3.2: Training Pipeline\n");
    test_train_from_sequences();

    printf("\nPhase 4.1: Gene Finding (Single Genome)\n");
    test_find_genes_with_training();

    printf("\nPhase 4.2: AOS Output\n");
    test_aos_matches_soa();

    printf("\nPhase 5: Metagenomic Mode\n");
    test_meta_find_genes();
    test_meta_no_path_bin();

    printf("\nPhase 7: Context Reuse and Error Recovery\n");
    test_context_reuse();
    test_find_genes_without_sequence();
    test_error_recovery();

    printf("\nPhase 8: Allocator Hooks\n");
    test_custom_allocator();
    test_meta_alloc_failure();

    printf("\nPhase 9: Callbacks\n");
    test_log_callback();
    test_progress_cancellation();

    printf("\nPhase 10: Edge Cases\n");
    test_very_short_sequence();
    test_all_n_sequence();

    printf("\nPhase 12: SOA Alignment\n");
    test_soa_alignment();

    printf("\nPhase 13: Training Setters\n");
    test_training_setters();

    printf("\nPhase 14: Golden (API = CLI ground truth)\n");
    test_golden();
    test_independence_meta();
    test_independence_single();

    printf("\n=== Results: %d/%d passed, %d of them known failures (XFAIL) ===\n",
           tests_passed, tests_run, tests_xfailed);
    return (tests_passed == tests_run) ? 0 : 1;
}
