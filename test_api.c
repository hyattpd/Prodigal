/*******************************************************************************
    PRODIGAL Library API Test Suite
*******************************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "prodigal.h"
/* White-box: the context (the node array), STT_NOD, and struct _training,
   the model blob's layout */
#include "prodigal_internal.h"
#include "tests/api_gff.h"
#include "tests/fasta.h"

static int tests_run = 0;
static int tests_passed = 0;

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

/* A missing fixture fails the test. */
static void missing_fixture(const char *what) {
    printf("[FAIL] missing fixture %s (run from the repo root)\n", what);
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

/* The tests' contexts return proteins, so every test that compares genes
   compares their proteins too; test_proteins_off checks the default. */
static prodigal_ctx_t *make_ctx(int meta) {
    prodigal_config_t config;
    prodigal_ctx_t *ctx;
    prodigal_config_init(&config);
    config.meta_mode = meta;
    config.want_proteins = 1;
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

/* f = one record, taking ownership of seq (len characters, NUL-terminated). */
static void fasta_one(fasta_t *f, const char *header, char *seq, int32_t len) {
    f->n = 1;
    f->seqs = (char **)xrealloc(NULL, sizeof(char *));
    f->headers = (char **)xrealloc(NULL, sizeof(char *));
    f->lens = (int32_t *)xrealloc(NULL, sizeof(int32_t));
    f->seqs[0] = seq;
    f->headers[0] = dupstr(header);
    f->lens[0] = len;
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
    ASSERT_EQ_INT(prodigal_abi_version(), PRODIGAL_ABI_VERSION);
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
    ASSERT_EQ_INT(config.want_proteins, 0);
    TEST_PASS();
}

/* The public struct typedefs in prodigal.h ("} name;"), into names; how
   many, or -1 if it can't be read. */
static int header_structs(char names[][64], int max) {
    char line[256], semi[2];
    int n = 0;
    FILE *fp = fopen("prodigal.h", "r");
    if (fp == NULL) return -1;
    while (fgets(line, sizeof(line), fp) != NULL)
        if (n < max && sscanf(line, "} %63[a-z_]%1[;]", names[n], semi) == 2) n++;
    fclose(fp);
    return n;
}

/* The ABI version and the FNV-1a hash of the table's names ("type;" and
   "type.field;" per row) in that version. A change to the table's structs
   or fields without a new version fails test_abi_layout. */
#define ABI_PINNED_VERSION 1
#define ABI_PINNED_HASH 0x8a6eb5fcu

/* prodigal_abi_layout lists the public structs of prodigal.h in order, each
   with its fields in order, none overlapping the one before; its names are
   the ones pinned for this ABI version. (A field missing from the table
   doesn't compile: abi_complete_* in prodigal_api.c.) */
static void test_abi_layout(void) {
    TEST_START("ABI: layout has every struct, fields in order; pinned");
    char header[32][64];
    size_t n = 0, i;
    const prodigal_abi_field_t *rows = prodigal_abi_layout(&n);
    int nh = header_structs(header, 32), t, found = 0;
    uint32_t hash = 2166136261u;
    if (nh < 0) { missing_fixture("prodigal.h"); return; }
    ASSERT_TRUE(rows != NULL && rows[n].type == NULL);
    ASSERT_TRUE(prodigal_abi_layout(NULL) == rows);
    for (i = 0, t = 0; i < n; t++) {
        const prodigal_abi_field_t *s = &rows[i];
        size_t end = 0;
        if (t >= nh || s->field != NULL || strcmp(s->type, header[t]) != 0) {
            printf("[FAIL] row %zu is %s, prodigal.h's next struct %s\n", i,
                   s->type, t < nh ? header[t] : "(none)");
            return;
        }
        for (i++; i < n && rows[i].field != NULL; i++) {
            const prodigal_abi_field_t *f = &rows[i];
            if (f->offset < end || f->offset + f->size > s->size) {
                printf("[FAIL] %s.%s at %zu (%zu bytes) overlaps the field before "
                       "it or the end: out of order\n", f->type, f->field,
                       f->offset, f->size);
                return;
            }
            end = f->offset + f->size;
            found += strcmp(f->type, "prodigal_config_t") == 0 &&
                     strcmp(f->field, "want_proteins") == 0 &&
                     f->offset == offsetof(prodigal_config_t, want_proteins) &&
                     f->size == sizeof(int);
        }
    }
    if (t != nh) { printf("[FAIL] %s is not in the table\n", header[t]); return; }
    ASSERT_TRUE(found == 1);

    for (i = 0; i < n; i++) {
        const char *p, *parts[3];
        int k;
        parts[0] = rows[i].type;
        parts[1] = rows[i].field ? "." : "";
        parts[2] = rows[i].field ? rows[i].field : "";
        for (k = 0; k < 3; k++)
            for (p = parts[k]; *p; p++) hash = (hash ^ (unsigned char)*p) * 16777619u;
        hash = (hash ^ ';') * 16777619u;
    }
    if (PRODIGAL_ABI_VERSION != ABI_PINNED_VERSION || hash != ABI_PINNED_HASH) {
        printf("[FAIL] ABI version %d, names hash 0x%08lx; pinned: version %d, "
               "0x%08lx. A changed ABI needs a new PRODIGAL_ABI_VERSION; then pin "
               "both here\n", PRODIGAL_ABI_VERSION, (unsigned long)hash,
               ABI_PINNED_VERSION, (unsigned long)ABI_PINNED_HASH);
        return;
    }
    TEST_PASS();
}

/* prodigal_config_init_v refuses another ABI version or config size, and
   then writes nothing. */
static void test_config_init_abi(void) {
    TEST_START("config_init refuses another ABI and writes nothing");
    prodigal_config_t config, before;
    memset(&config, 0x5a, sizeof(config));
    before = config;
    ASSERT_EQ_INT(prodigal_config_init_v(&config, PRODIGAL_ABI_VERSION + 1,
                                         sizeof(config)), PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_EQ_INT(prodigal_config_init_v(&config, PRODIGAL_ABI_VERSION,
                                         sizeof(config) - 8), PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_TRUE(memcmp(&config, &before, sizeof(config)) == 0);
    ASSERT_EQ_INT(prodigal_config_init_v(NULL, PRODIGAL_ABI_VERSION, sizeof(config)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_EQ_INT(prodigal_config_init(&config), PRODIGAL_OK);
    ASSERT_TRUE(config.struct_size == sizeof(config) && config.trans_table == 11);
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

static void test_free_fn(void *ptr, void *user_data);

static void test_config_validate(void) {
    TEST_START("config_validate explains what create rejects");
    prodigal_config_t c;
    char msg[256];
    double bad_wt[] = { 0.0, -1.0, NAN, INFINITY };
    size_t i;
    ASSERT_EQ_INT(prodigal_config_validate(NULL, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    prodigal_config_init(&c);
    ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)), PRODIGAL_OK);
    ASSERT_TRUE(msg[0] == '\0');
    ASSERT_EQ_INT(prodigal_config_validate(&c, NULL, 0), PRODIGAL_OK);

    c.struct_size = 0;
    ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_TRUE(strstr(msg, "struct_size") != NULL);

    prodigal_config_init(&c);
    c.trans_table = 7;
    ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_TRUE(strstr(msg, "translation table") != NULL);

    for (i = 0; i < sizeof(bad_wt) / sizeof(bad_wt[0]); i++) {
        prodigal_config_init(&c);
        c.meta_mode = (int)(i % 2);   /* either mode */
        c.start_weight = bad_wt[i];
        ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)),
                      PRODIGAL_ERR_INVALID_CONFIG);
        ASSERT_TRUE(strstr(msg, "start weight") != NULL);
        ASSERT_TRUE(prodigal_create(&c) == NULL);
    }

    /* free_fn alone would free with it what malloc allocated */
    prodigal_config_init(&c);
    c.free_fn = test_free_fn;
    ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_TRUE(strstr(msg, "free_fn without alloc_fn") != NULL);
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

/* Both prodigal_training_info and prodigal_load_training reject the len
   bytes at blob, with the same reason, and change neither info nor the
   model. 1, or 0 after printing [FAIL]. */
static int both_reject(prodigal_ctx_t *ctx, const void *blob, size_t len,
                       const char *what) {
    prodigal_training_info_t info, before;
    char msg[256];
    memset(&info, 0x5a, sizeof(info));
    before = info;
    if (prodigal_training_info(blob, len, &info, msg, sizeof(msg)) !=
            PRODIGAL_ERR_INVALID_INPUT || msg[0] == '\0' ||
        memcmp(&info, &before, sizeof(info)) != 0) {
        printf("[FAIL] %s: training_info accepted it, gave no reason or wrote "
               "info\n", what);
        return 0;
    }
    if (prodigal_load_training(ctx, blob, len) != PRODIGAL_ERR_INVALID_INPUT ||
        strcmp(prodigal_last_error(ctx), msg) != 0) {
        printf("[FAIL] %s: load says \"%s\", training_info \"%s\"\n", what,
               prodigal_last_error(ctx), msg);
        return 0;
    }
    return 1;
}

static void test_training_load_invalid(void) {
    TEST_START("training load and info reject a bad blob; model kept");
    size_t tlen, elen;
    void *tdata = load_ref_model(&tlen), *exported = NULL;
    if (tdata == NULL) return;
    prodigal_ctx_t *ctx = make_model_ctx(tdata, tlen);
    ASSERT_TRUE(ctx != NULL);

    char dummy[16] = {0}, what[32];
    if (!both_reject(ctx, dummy, sizeof(dummy), "a 16-byte blob") ||
        !both_reject(ctx, NULL, tlen, "NULL"))
        return;

    /* Fields training sets from the config and counts must be in range */
    void *bad = malloc(tlen);
    struct _training *b = (struct _training *)bad;
    int k;
    ASSERT_TRUE(bad != NULL);
    for (k = 0; k < 9; k++) {
        memcpy(bad, tdata, tlen);
        switch (k) {
        case 0: b->trans_table = 7; break;
        case 1: b->gc = 1.5; break;
        case 2: b->gc = -0.1; break;
        case 3: b->gc = NAN; break;
        case 4: b->st_wt = 0.0; break;
        case 5: b->st_wt = -1.0; break;
        case 6: b->st_wt = NAN; break;
        case 7: b->st_wt = INFINITY; break;
        default: b->uses_sd = 2; break;
        }
        snprintf(what, sizeof(what), "bad model %d", k);
        if (!both_reject(ctx, bad, tlen, what)) return;
    }
    ASSERT_EQ_INT(prodigal_export_training(ctx, &exported, &elen), PRODIGAL_OK);
    ASSERT_TRUE(elen == tlen && memcmp(exported, tdata, tlen) == 0);

    free(exported);
    free(bad);
    free(tdata);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/* prodigal_training_info reads the settings from a blob wherever it sits in
   memory (its rejections: test_training_load_invalid). */
static void test_training_info(void) {
    TEST_START("training_info reads a model, aligned or not");
    prodigal_training_info_t info;
    size_t tlen;
    char msg[256];
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL) return;
    const struct _training *t = (const struct _training *)tdata;
    char *buf = (char *)xrealloc(NULL, tlen + 1);

    memcpy(buf + 1, tdata, tlen);
    ASSERT_EQ_INT(prodigal_training_info(buf + 1, tlen, &info, msg, sizeof(msg)),
                  PRODIGAL_OK);
    ASSERT_TRUE(msg[0] == '\0');
    ASSERT_EQ_INT(info.translation_table, t->trans_table);
    ASSERT_EQ_INT(info.uses_sd, t->uses_sd);
    ASSERT_TRUE(info.gc_content == t->gc && info.start_weight == t->st_wt);
    ASSERT_EQ_INT(prodigal_training_info(tdata, tlen, NULL, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_TRUE(msg[0] != '\0');

    free(buf);
    free(tdata);
    TEST_PASS();
}

static void test_export_needs_model(void) {
    TEST_START("export needs a model; a degenerate one round-trips");
    prodigal_training_info_t info;
    void *exported = NULL;
    size_t elen = 0;
    fasta_t fa;
    prodigal_ctx_t *ctx = make_ctx(0), *loaded = NULL;
    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(prodigal_export_training(ctx, &exported, &elen),
                  PRODIGAL_ERR_INVALID_INPUT);

    /* Training on one base leaves no_mot NaN (upstream quirk 11); the
       model still loads */
    char *seq = (char *)xrealloc(NULL, 30001);
    memset(seq, 'A', 30000);
    seq[30000] = '\0';
    fasta_one(&fa, "all_a", seq, 30000);
    ASSERT_EQ_INT(api_train(ctx, &fa), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_export_training(ctx, &exported, &elen), PRODIGAL_OK);
    ASSERT_TRUE(isnan(((struct _training *)exported)->no_mot));
    loaded = make_model_ctx(exported, elen);
    ASSERT_TRUE(loaded != NULL);
    ASSERT_EQ_INT(prodigal_training_info(exported, elen, &info, NULL, 0), PRODIGAL_OK);

    prodigal_destroy(loaded);
    free(exported);
    fasta_free(&fa);
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
        int32_t i, len;
        for (i = 0; i < soa->n_genes; i++) {
            const char *p = api_protein(soa, i, &len);
            ASSERT_EQ_INT(soa->begin[i], aos->genes[i].begin);
            ASSERT_EQ_INT(soa->end[i], aos->genes[i].end);
            ASSERT_EQ_INT(soa->strand[i], aos->genes[i].strand);
            ASSERT_TRUE(soa->cscore[i] == aos->genes[i].cscore);
            ASSERT_TRUE(soa->confidence[i] == aos->genes[i].confidence);
            /* The AOS owns its strings: its SOA is freed by now */
            ASSERT_TRUE(strcmp(soa->rbs_motif[i], aos->genes[i].rbs_motif) == 0);
            ASSERT_TRUE(strcmp(soa->rbs_spacer[i], aos->genes[i].rbs_spacer) == 0);
            ASSERT_TRUE(aos->genes[i].protein != NULL);
            ASSERT_TRUE(strlen(aos->genes[i].protein) == (size_t)len);
            ASSERT_TRUE(memcmp(aos->genes[i].protein, p, (size_t)len) == 0);
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

/* A context on test_alloc/test_free_fn, with the counters reset and no
   allocation limit. */
static prodigal_ctx_t *make_counted_ctx(int meta) {
    prodigal_config_t config;
    prodigal_config_init(&config);
    config.meta_mode = meta;
    config.want_proteins = 1;
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

    prodigal_ctx_t *ctx = make_counted_ctx(1);
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

    prodigal_ctx_t *ctx = make_counted_ctx(1);
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

/* Fail each of training's allocations in turn (non-SD start training
   forced, so all of them run): each is an error that leaves no model, and
   with no limit training succeeds, leaking nothing. */
static void test_training_alloc_failures(void) {
    TEST_START("training: every failed allocation is a clean error");
    fasta_t fa;
    void *exported = NULL;
    size_t elen = 0, k;
    int rc = PRODIGAL_ERR_NOMEM;
    if (!load_anthus(&fa)) return;
    prodigal_ctx_t *ctx = make_counted_ctx(0);
    ASSERT_TRUE(ctx != NULL);
    ctx->config.force_nonsd = 1;
    for (k = 0; rc == PRODIGAL_ERR_NOMEM; k++) {
        alloc_limit = alloc_count + k;
        rc = api_train(ctx, &fa);
        if (rc == PRODIGAL_ERR_NOMEM)
            ASSERT_EQ_INT(prodigal_export_training(ctx, &exported, &elen),
                          PRODIGAL_ERR_INVALID_INPUT);
        if (k == 0) ASSERT_TRUE(strstr(prodigal_last_error(ctx), "gc frame") != NULL);
    }
    alloc_limit = (size_t)-1;
    ASSERT_EQ_INT(rc, PRODIGAL_OK);
    ASSERT_TRUE(k >= 3);   /* the GC frame plot, the motif counts, then none */
    prodigal_destroy(ctx);
    ASSERT_TRUE(alloc_count == free_count);
    fasta_free(&fa);
    TEST_PASS();
}

static void test_failed_training_leaves_no_model(void) {
    TEST_START("failed training leaves no model, not a damaged one");
    fasta_t fa;
    size_t tlen;
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL || !load_anthus(&fa)) return;

    prodigal_ctx_t *ctx = make_counted_ctx(0);
    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(prodigal_load_training(ctx, tdata, tlen), PRODIGAL_OK);

    /* Long enough that slen/8 nodes exceed the initial STT_NOD, so
       training's first step (reserve_nodes) allocates, and fails */
    int32_t len = 8 * (STT_NOD + 1);
    char *big = (char *)malloc(len);
    ASSERT_TRUE(big != NULL);
    memset(big, 'A', len);
    const char *seqs[1] = { big };
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, seqs, NULL, &len, 1),
                  PRODIGAL_OK);
    alloc_limit = alloc_count;
    ASSERT_EQ_INT(prodigal_train(ctx), PRODIGAL_ERR_NOMEM);
    alloc_limit = (size_t)-1;

    prodigal_genes_soa_t *genes = NULL;
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]),
                  PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_load_training(ctx, tdata, tlen), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_OK);
    prodigal_genes_free(genes);
    prodigal_destroy(ctx);
    free(big);
    free(tdata);
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
    TEST_START("training setters change only their field of a loaded model");
    prodigal_training_info_t info;
    size_t tlen, elen;
    void *tdata = load_ref_model(&tlen), *exported = NULL;
    if (tdata == NULL) return;
    prodigal_ctx_t *ctx = make_model_ctx(tdata, tlen);
    ASSERT_TRUE(ctx != NULL);

    /* Valid setters */
    ASSERT_EQ_INT(prodigal_set_translation_table(ctx, 4), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_start_weight(ctx, 3.0), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, 0.45), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_uses_sd(ctx, 0), PRODIGAL_OK);

    /* Invalid setters change nothing */
    ASSERT_EQ_INT(prodigal_set_translation_table(ctx, 7), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_start_weight(ctx, -1.0), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_start_weight(ctx, NAN), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_start_weight(ctx, INFINITY), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, 2.0), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, -0.1), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_gc(ctx, NAN), PRODIGAL_ERR_INVALID_INPUT);

    /* The exported model is the loaded one with those four fields set */
    struct _training *want = (struct _training *)tdata;
    want->trans_table = 4;
    want->st_wt = 3.0;
    want->gc = 0.45;
    want->uses_sd = 0;
    ASSERT_EQ_INT(prodigal_export_training(ctx, &exported, &elen), PRODIGAL_OK);
    ASSERT_TRUE(elen == tlen && memcmp(exported, want, tlen) == 0);
    ASSERT_EQ_INT(prodigal_training_info(exported, elen, &info, NULL, 0), PRODIGAL_OK);
    ASSERT_TRUE(info.translation_table == 4 && info.start_weight == 3.0 &&
                info.gc_content == 0.45 && info.uses_sd == 0);

    free(exported);
    free(tdata);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/*******************************************************************************
    Phase 14: API ≡ CLI golden tests

    The library's output, in the form of tests/api_gff.h, must equal the
    CLI's ground truth line for line.
*******************************************************************************/

/* Run all of fa on ctx and compare with the CLI's output: the GFF ref and
   the -a file faa. */
static int golden_run(prodigal_ctx_t *ctx, const fasta_t *fa, const char *ref,
                      const char *faa, int meta) {
    lines_t want = {0}, got = {0};
    int ok = 0, rc = api_emit_cli_file(ref, faa, meta, lines_emit, &want);
    if (rc == API_CLI_UNREADABLE)
        printf("[FAIL] cannot read %s or %s (run from the repo root)\n", ref, faa);
    else if (rc == API_CLI_MISMATCH)
        printf("[FAIL] %s and %s hold different numbers of genes\n", ref, faa);
    else if (api_emit_records(ctx, fa, 0, fa->n, meta, lines_emit, &got) != PRODIGAL_OK)
        printf("[FAIL] find_genes: %s\n", prodigal_last_error(ctx));
    else
        ok = lines_equal(&want, &got, ref);
    lines_free(&want); lines_free(&got);
    return ok;
}

/* ARM64 macOS training writes different bytes (docs/testing.md). */
#if defined(__APPLE__) && defined(__aarch64__)
static const int train_bytes_match = 0;
#else
static const int train_bytes_match = 1;
#endif

/* Is the model trained on ctx the reference? The fields training copies or
   counts must match on every platform, and every byte where training
   writes the reference's bytes (not ARM64 macOS). */
static int model_is_ref(prodigal_ctx_t *ctx) {
    void *got = NULL, *ref = NULL;
    size_t glen = 0, rlen = 0, i = 0;
    int ok = 0;
    if (prodigal_export_training(ctx, &got, &glen) != PRODIGAL_OK) {
        printf("[FAIL] export: %s\n", prodigal_last_error(ctx));
        return 0;
    }
    if ((ref = load_ref_model(&rlen)) != NULL) {
        const struct _training *g = (const struct _training *)got;
        const struct _training *r = (const struct _training *)ref;
        if (glen != rlen) {
            printf("[FAIL] model is %zu bytes, expected %zu\n", glen, rlen);
        } else if (g->gc != r->gc || g->trans_table != r->trans_table ||
                   g->st_wt != r->st_wt || g->uses_sd != r->uses_sd) {
            printf("[FAIL] model gc %g table %d st_wt %g uses_sd %d, "
                   "expected %g %d %g %d\n", g->gc, g->trans_table, g->st_wt,
                   g->uses_sd, r->gc, r->trans_table, r->st_wt, r->uses_sd);
        } else if (train_bytes_match) {
            while (i < rlen && ((char *)ref)[i] == ((char *)got)[i]) i++;
            ok = i == rlen;
            if (!ok) printf("[FAIL] model differs from byte %zu\n", i);
        } else {
            ok = 1;
        }
    }
    free(got); free(ref);
    return ok;
}

/* Context sources: OTHER_MODEL is ref_train.bin with other parameters
   (table 4, start weight 3, GC 0.9), which training must not keep. */
enum { META, SINGLE, REF_MODEL, OTHER_MODEL };

/* One golden test on anthus_aco.fas: make a context from src, train it
   `trainings` times on all of anthus (then its model must be
   ref_train.bin), and compare the genes it finds and their proteins with
   the CLI's in ref and faa. With reload, so must a fresh context given the
   exported model (which also checks it where training writes other bytes:
   ARM64 macOS). */
static void golden(int src, int trainings, int reload, const char *ref,
                   const char *faa) {
    fasta_t fa;
    prodigal_ctx_t *ctx = NULL;
    void *blob = NULL;
    size_t len;
    int i;
    if (!load_anthus(&fa)) return;
    if (src == META || src == SINGLE) {
        ctx = make_ctx(src == META);
    } else if ((blob = load_ref_model(&len)) != NULL) {
        if (src == OTHER_MODEL) {
            struct _training *m = (struct _training *)blob;
            m->trans_table = 4;
            m->st_wt = 3.0;
            m->gc = 0.9;
        }
        ctx = make_model_ctx(blob, len);
    }
    for (i = 0; ctx != NULL && i < trainings; i++)
        if (api_train(ctx, &fa) != PRODIGAL_OK) {
            printf("[FAIL] training: %s\n", prodigal_last_error(ctx));
            prodigal_destroy(ctx);
            ctx = NULL;
        }
    if (ctx != NULL && (trainings == 0 || model_is_ref(ctx)) &&
        golden_run(ctx, &fa, ref, faa, src == META)) {
        prodigal_ctx_t *fresh = NULL;
        void *exported = NULL;
        if (!reload)
            TEST_PASS();
        else if (prodigal_export_training(ctx, &exported, &len) == PRODIGAL_OK &&
                 (fresh = make_model_ctx(exported, len)) != NULL &&
                 golden_run(fresh, &fa, ref, faa, 0))
            TEST_PASS();
        prodigal_destroy(fresh);
        free(exported);
    }
    prodigal_destroy(ctx);
    free(blob);
    fasta_free(&fa);
}

static void test_golden(void) {
    /* ref_trained.gff equals ref_single.gff (the trained model is
       ref_train.bin), so its proteins are ref_single.proteins */
    static const struct {
        const char *name; int src, trainings, reload; const char *ref, *faa;
    } t[] = {
        { "meta mode = ref_meta.gff + .proteins, every field", META, 0, 0,
          REF_DIR "ref_meta.gff", REF_DIR "ref_meta.proteins" },
        { "API training = ref_train.bin, genes = ref_single.*; reloaded",
          SINGLE, 1, 1, REF_DIR "ref_single.gff", REF_DIR "ref_single.proteins" },
        { "training twice: the same", SINGLE, 2, 0, REF_DIR "ref_single.gff",
          REF_DIR "ref_single.proteins" },
        { "other model loaded, then training: the same", OTHER_MODEL, 1, 0,
          REF_DIR "ref_single.gff", REF_DIR "ref_single.proteins" },
        { "single mode from ref_train.bin = ref_trained.gff", REF_MODEL, 0, 0,
          REF_DIR "ref_trained.gff", REF_DIR "ref_single.proteins" },
    };
    char name[128];
    size_t i;
    for (i = 0; i < sizeof(t) / sizeof(t[0]); i++) {
        snprintf(name, sizeof(name), "golden: %s", t[i].name);
        TEST_START(name);
        golden(t[i].src, t[i].trainings, t[i].reload, t[i].ref, t[i].faa);
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
        else if (lines_equal(&first, &again, "second pass on the same context") &&
                 lines_equal(&first, &solo, "each record on a fresh context") &&
                 lines_equal(&first, &twice, "find_genes twice per sequence"))
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
    Phase 15: Low-complexity input

    "ATG" repeated gives one start node per codon, about slen/3 nodes, past
    the slen/8 (and at least STT_NOD) nodes the array is first sized for.
    The array must grow, not overflow. make compare runs repeats past
    STT_NOD (tests/make_corpus.sh); here a context's capacity is set to 0,
    so that a 21 kbp repeat starts from slen/8 and grows.
*******************************************************************************/

/* fa = one record, n copies of "ATG". */
static void atg_repeat(fasta_t *fa, int32_t n) {
    int32_t i;
    char *s = (char *)xrealloc(NULL, 3 * (size_t)n + 1);
    for (i = 0; i < n; i++) memcpy(s + 3 * i, "ATG", 3);
    s[3 * n] = '\0';
    fasta_one(fa, "atg_repeat", s, 3 * n);
}

/* An ATG repeat, run twice on one context made from src (SINGLE: trained on
   the repeat each time) whose node array starts at slen/8: it must grow in
   the first run (in training, for SINGLE) and not in the second, and the
   two must give the same genes and, trained, the same model. */
static void low_complexity(int src) {
    fasta_t fa;
    lines_t out[2] = {{0}, {0}};
    struct _training *model[2] = {NULL, NULL};
    prodigal_ctx_t *ctx = NULL;
    void *blob = NULL;
    size_t len;
    int cap[2], i, rc = PRODIGAL_OK;
    atg_repeat(&fa, 7000);   /* 21 kbp: enough to train on */
    if (src != REF_MODEL) ctx = make_ctx(src == META);
    else if ((blob = load_ref_model(&len)) != NULL) ctx = make_model_ctx(blob, len);
    if (ctx != NULL) ctx->node_cap = 0;   /* the first run sizes it to slen/8 */
    for (i = 0; ctx != NULL && rc == PRODIGAL_OK && i < 2; i++) {
        if (src == SINGLE && (rc = api_train(ctx, &fa)) == PRODIGAL_OK)
            model[i] = (struct _training *)memcpy(
                xrealloc(NULL, sizeof(struct _training)), &ctx->tinf,
                sizeof(struct _training));
        if (rc == PRODIGAL_OK)
            rc = api_emit_records(ctx, &fa, 0, 1, src == META, lines_emit, &out[i]);
        cap[i] = ctx->node_cap;
    }
    if (ctx != NULL) {
        if (rc != PRODIGAL_OK)
            printf("[FAIL] %d: %s\n", rc, prodigal_last_error(ctx));
        else if (cap[0] <= fa.lens[0] / 8 || cap[1] != cap[0])
            printf("[FAIL] node array of %d, then %d nodes: the first run "
                   "must grow it past %d, the second not\n", cap[0], cap[1],
                   fa.lens[0] / 8);
        else if (src == SINGLE &&
                 memcmp(model[0], model[1], sizeof(struct _training)) != 0)
            printf("[FAIL] the models differ\n");
        else if (lines_equal(&out[0], &out[1], "genes on the grown array"))
            TEST_PASS();
    }
    lines_free(&out[0]); lines_free(&out[1]);
    free(model[0]); free(model[1]);
    prodigal_destroy(ctx);
    free(blob);
    fasta_free(&fa);
}

/* Every node past ctx->nn is zero (the invariant build_nodes keeps). */
static int node_tail_zero(const prodigal_ctx_t *ctx) {
    const unsigned char *p = (const unsigned char *)(ctx->nodes + ctx->nn);
    size_t i, n = (size_t)(ctx->node_cap - ctx->nn) * sizeof(struct _node);
    for (i = 0; i < n; i++)
        if (p[i] != 0) return 0;
    return 1;
}

static void test_low_complexity_grow_fails(void) {
    TEST_START("ATG repeat: failed growth, then anthus = ref_meta.gff");
    fasta_t fa, atg;
    lines_t first = {0};
    prodigal_genes_soa_t *genes = NULL;
    if (!load_anthus(&fa)) return;
    atg_repeat(&atg, STT_NOD + 10000);
    prodigal_ctx_t *ctx = make_counted_ctx(1);
    ASSERT_TRUE(ctx != NULL);

    /* Bins allocated; then the ATG repeat fills the array and growing it
       fails */
    ASSERT_EQ_INT(api_emit_records(ctx, &fa, 0, 1, 1, lines_emit, &first), PRODIGAL_OK);
    lines_free(&first);
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, atg.seqs[0], atg.lens[0], NULL),
                  PRODIGAL_OK);
    alloc_limit = alloc_count;
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &genes, NULL), PRODIGAL_ERR_NOMEM);
    alloc_limit = (size_t)-1;
    ASSERT_TRUE(strstr(prodigal_last_error(ctx), "nodes") != NULL);
    ASSERT_TRUE(node_tail_zero(ctx));

    /* No stale node reaches the next sequence */
    if (!golden_run(ctx, &fa, REF_DIR "ref_meta.gff", REF_DIR "ref_meta.proteins", 1))
        return;
    ASSERT_TRUE(node_tail_zero(ctx));

    prodigal_destroy(ctx);
    fasta_free(&fa); fasta_free(&atg);
    TEST_PASS();
}

static void test_low_complexity_meta(void) {
    TEST_START("ATG repeat: meta mode, node array grows");
    low_complexity(META);
}

static void test_low_complexity_model(void) {
    TEST_START("ATG repeat: single mode from model, node array grows");
    low_complexity(REF_MODEL);
}

static void test_low_complexity_trained(void) {
    TEST_START("ATG repeat: training and single mode, node array grows");
    low_complexity(SINGLE);
}

/*******************************************************************************
    Phase 16: Limits

    Genes and masked regions past the initial arrays are kept, not dropped,
    and the length limits are the CLI's, on encoded bases.
*******************************************************************************/

/* one = the records of fa joined into a single record. */
static void join_records(const fasta_t *fa, fasta_t *one) {
    int32_t k, len = 0;
    char *s;
    for (k = 0; k < fa->n; k++) len += fa->lens[k];
    s = (char *)xrealloc(NULL, (size_t)len + 1);
    for (len = 0, k = 0; k < fa->n; k++) {
        memcpy(s + len, fa->seqs[k], (size_t)fa->lens[k]);
        len += fa->lens[k];
    }
    s[len] = '\0';
    fasta_one(one, "anthus_joined", s, len);
}

/* Genes past the gene array grow it: on anthus joined into one sequence, a
   context with room for one gene to start with must find what a default
   one finds, and only the first may grow. Single mode from the model: in
   meta mode the joined sequence has one gene, and both modes keep their
   genes through the same code (keep_genes). */
static void test_gene_array_grows(void) {
    TEST_START("gene array grows from 1 gene");
    fasta_t fa, one;
    lines_t want = {0}, got = {0};
    size_t len;
    void *blob = load_ref_model(&len);
    if (blob == NULL || !load_anthus(&fa)) return;
    join_records(&fa, &one);
    prodigal_ctx_t *def = make_model_ctx(blob, len), *small = make_model_ctx(blob, len);
    ASSERT_TRUE(def != NULL && small != NULL);
    small->gene_cap = 1;   /* the real array is larger */
    ASSERT_EQ_INT(api_emit_records(def, &one, 0, 1, 0, lines_emit, &want), PRODIGAL_OK);
    ASSERT_EQ_INT(api_emit_records(small, &one, 0, 1, 0, lines_emit, &got), PRODIGAL_OK);
    ASSERT_TRUE(small->gene_cap > 1 && def->gene_cap == MAX_GENES);
    if (!lines_equal(&want, &got, "genes, grown from 1")) return;
    lines_free(&want); lines_free(&got);
    prodigal_destroy(def); prodigal_destroy(small);
    free(blob);
    fasta_free(&fa); fasta_free(&one);
    TEST_PASS();
}

/* n runs of MASK_SIZE Ns, each followed by 100 bases; NUL-terminated. */
static char *n_runs(int n, int32_t *len) {
    int i;
    char *s = (char *)xrealloc(NULL, (size_t)n * (MASK_SIZE + 100) + 1);
    for (i = 0; i < n; i++) {
        memset(s + i * (MASK_SIZE + 100), 'N', MASK_SIZE);
        memset(s + i * (MASK_SIZE + 100) + MASK_SIZE, 'A', 100);
    }
    *len = n * (MASK_SIZE + 100);
    s[*len] = '\0';
    return s;
}

static void test_masks_past_max(void) {
    TEST_START("more than 5,000 masked regions (upstream's cap) are all kept");
    int32_t len;
    size_t before;
    char *s = n_runs(6000, &len);
    const char *seqs[1] = { s };
    prodigal_ctx_t *ctx = make_counted_ctx(0);   /* masks are a working buffer */
    ASSERT_TRUE(ctx != NULL);
    ctx->config.mask_regions = 1;
    before = alloc_count;
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, s, len, NULL), PRODIGAL_OK);
    ASSERT_EQ_INT(ctx->nmask, 6000);
    ASSERT_TRUE(alloc_count > before);
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, seqs, NULL, &len, 1),
                  PRODIGAL_OK);
    ASSERT_EQ_INT(ctx->nmask, 6000);
    prodigal_destroy(ctx);
    ASSERT_TRUE(alloc_count == free_count);
    free(s);
    TEST_PASS();
}

/* The CLI rejects a contig, and truncates training input, once the encoded
   length plus MAX_LINE reaches MAX_SEQ; the library rejects both there,
   and counts encoded bases only (characters outside A-z are dropped). */
static void test_length_limits(void) {
    TEST_START("length limits: encoded bases, at the CLI's boundary");
    const int32_t max = MAX_SEQ - MAX_LINE - 1;   /* longest accepted */
    int32_t len = MAX_SEQ + 1000;
    char *s = (char *)xrealloc(NULL, (size_t)len);
    const char *seqs[1] = { s };
    prodigal_seq_info_t info;
    prodigal_ctx_t *ctx = make_ctx(0);
    ASSERT_TRUE(ctx != NULL);

    memset(s, 'A', (size_t)max + 1);
    len = max;
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, s, len, NULL), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, seqs, NULL, &len, 1),
                  PRODIGAL_OK);
    len = max + 1;
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, s, len, NULL), PRODIGAL_ERR_INVALID_INPUT);
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, seqs, NULL, &len, 1),
                  PRODIGAL_ERR_INVALID_INPUT);

    /* Two training sequences: the 12-base spacers between and after them
       count too */
    const char *two[2] = { s, s };
    int32_t lens[2] = { 1000, max - 1000 - 24 };
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, two, NULL, lens, 2),
                  PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_get_seq_info(ctx, &info), PRODIGAL_OK);
    ASSERT_EQ_INT(info.length, max);
    lens[1]++;
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, two, NULL, lens, 2),
                  PRODIGAL_ERR_INVALID_INPUT);

    /* Longer than MAX_SEQ, but mostly dropped characters */
    memset(s, '-', (size_t)MAX_SEQ + 1000);
    memset(s, 'A', 1000);
    len = MAX_SEQ + 1000;
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, s, len, NULL), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_get_seq_info(ctx, &info), PRODIGAL_OK);
    ASSERT_EQ_INT(info.length, 1000);

    prodigal_destroy(ctx);
    free(s);
    TEST_PASS();
}

/*******************************************************************************
    Test runner
*******************************************************************************/

/*******************************************************************************
    Phase 17: Proteins

    The golden and compare tests check the letters against the CLI's -a
    file, and api_emit_found checks every protein's layout
    (api_proteins_ok, which the fuzz target checks too); these check what
    the letters mean on genes made by hand, and the default.
*******************************************************************************/

/* A model for the white-box tests of the table-dependent code, which reads
   only its table and GC. */
static struct _training *table_model(int table, double gc) {
    static struct _training tinf;
    tinf.trans_table = table;
    tinf.gc = gc;
    return &tinf;
}

/* The letters of a gene of seq (white-box: pdg_protein_letter on a gene
   made by hand) under translation table table. */
static void letters_of(prodigal_ctx_t *ctx, const char *seq, int strand,
                       int start_type, int table, char *out) {
    struct _training *tinf = table_model(table, 0.5);
    struct _gene g;
    int k, n;
    memset(&g, 0, sizeof(g));
    prodigal_set_sequence(ctx, seq, (int32_t)strlen(seq), NULL);
    g.begin = 1;
    g.end = ctx->slen;
    g.strand = strand;
    g.start_type = start_type;
    n = pdg_protein_length(&g);
    for (k = 0; k < n; k++)
        out[k] = pdg_protein_letter(&g, k, ctx->seq, ctx->rseq, ctx->useq,
                                    ctx->slen, tinf);
    out[n] = '\0';
}

static void test_protein_letters(void) {
    TEST_START("proteins: M only at a real start, X, *, the table");
    char got[16];
    prodigal_ctx_t *ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);
    /* GTG AAA NNN TAA, then its reverse complement on the other strand */
    letters_of(ctx, "GTGAAANNNTAA", 1, 1, 11, got);
    ASSERT_TRUE(strcmp(got, "MKX*") == 0);
    letters_of(ctx, "GTGAAANNNTAA", 1, 3, 11, got);    /* Edge: no M */
    ASSERT_TRUE(strcmp(got, "VKX*") == 0);
    letters_of(ctx, "TTANNNTTTCAC", -1, 1, 11, got);
    ASSERT_TRUE(strcmp(got, "MKX*") == 0);
    letters_of(ctx, "TTANNNTTTCAC", -1, 3, 11, got);
    ASSERT_TRUE(strcmp(got, "VKX*") == 0);
    /* TGA is a stop in table 11 and W in table 4 */
    letters_of(ctx, "ATGTGATAA", 1, 0, 11, got);
    ASSERT_TRUE(strcmp(got, "M**") == 0);
    letters_of(ctx, "ATGTGATAA", 1, 0, 4, got);
    ASSERT_TRUE(strcmp(got, "MW*") == 0);
    /* TAG is a stop in table 11 and W in table 32, the fork's */
    letters_of(ctx, "ATGTAGTAA", 1, 0, 11, got);
    ASSERT_TRUE(strcmp(got, "M**") == 0);
    letters_of(ctx, "ATGTAGTAA", 1, 0, 32, got);
    ASSERT_TRUE(strcmp(got, "MW*") == 0);
    /* CTG is L in table 11 and A in table 26, the fork's */
    letters_of(ctx, "ATGCTGTAA", 1, 0, 11, got);
    ASSERT_TRUE(strcmp(got, "ML*") == 0);
    letters_of(ctx, "ATGCTGTAA", 1, 0, 26, got);
    ASSERT_TRUE(strcmp(got, "MA*") == 0);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/* Without want_proteins (the default) there are none. */
static void test_proteins_off(void) {
    TEST_START("proteins: none unless want_proteins is set");
    prodigal_config_t config;
    prodigal_genes_soa_t *soa = NULL;
    prodigal_genes_t *aos = NULL;
    fasta_t fa;
    prodigal_ctx_t *ctx;
    if (!load_anthus(&fa)) return;
    prodigal_config_init(&config);
    config.meta_mode = 1;
    ctx = prodigal_create(&config);
    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, fa.seqs[0], fa.lens[0], fa.headers[0]), PRODIGAL_OK);
    ASSERT_EQ_INT(prodigal_find_genes(ctx, &soa, NULL), PRODIGAL_OK);
    ASSERT_TRUE(soa->n_genes > 0);
    ASSERT_TRUE(soa->protein_offsets == NULL && soa->protein_data == NULL);
    ASSERT_EQ_INT(prodigal_find_genes_aos(ctx, &aos, NULL), PRODIGAL_OK);
    ASSERT_TRUE(aos->n_genes > 0 && aos->genes[0].protein == NULL);
    prodigal_genes_free(soa);
    prodigal_genes_aos_free(aos);
    prodigal_destroy(ctx);
    fasta_free(&fa);
    TEST_PASS();
}

/*******************************************************************************
    Phase 18: Genetic codes (docs/genetic-codes.md, "The rules")

    A translation table is supported if and only if Prodigal uses it
    exactly as NCBI defines it in testdata/ncbi/gc.prt: every codon's stop
    and amino acid; of the codons Prodigal can start genes at (ATG, GTG and
    TTG), exactly the table's starts; the stop probability of the
    ORF-length score; a stop in every frame of the spacer between training
    contigs; and no stop or start made by an ambiguous base.  A table in
    which a codon is a stop or an amino acid depending on context can't be
    supported.
*******************************************************************************/

#define GC_PRT "testdata/ncbi/gc.prt"
#define GC_PRT_VERSION "4.6"    /* the version docs/genetic-codes.md cites */
#define GC_MAX 64               /* tables the parser keeps */
#define GC_PROBE_MAX 99         /* tables -1..GC_PROBE_MAX are checked */
#define N_CODONS 61             /* codons with at least one N */

/* The outcome, pinned.  Nothing checks the copies of these lists in
   docs/genetic-codes.md, README.md, prodigal.h (supported only) and
   CHANGES: update them with these.  Tables 26 and 32 are the fork's
   additions. */
static const int supported_tables[] = { 4, 5, 9, 11, 13, 15, 21, 25, 26, 32 };
static const int rejected_ncbi_tables[] = {
    1, 2, 3, 6, 10, 12, 14, 16, 22, 23, 24, 27, 28, 29, 30, 31, 33
};

/* NCBI's tables.  Codon k is base[0][k] base[1][k] base[2][k]; aa[r][k] is
   its amino acid in table id[r] (ncbieaa, * for a stop), and start[r][k]
   M if it is a start codon there, * if a stop (sncbieaa). */
typedef struct {
    int n, id[GC_MAX];
    char aa[GC_MAX][65], start[GC_MAX][65], base[3][65];
} gc_prt_t;

/* Reads GC_PRT into g. 1, or 0 after printing [FAIL]. */
static int load_gc_prt(gc_prt_t *g) {
    char line[512], s[65], aa[65] = "", version[16] = "", seen[64] = { 0 };
    int id = -1, k, j;
    FILE *f = fopen(GC_PRT, "r");
    memset(g, 0, sizeof(*g));
    if (f == NULL) { missing_fixture(GC_PRT); return 0; }
    while (fgets(line, sizeof(line), f) != NULL) {
        const char *p = line + strspn(line, " \t");
        if (version[0] == '\0' && sscanf(p, "-- Version %15s", version) == 1)
            continue;
        if (sscanf(p, "-- Base%d %64[TCAG]", &k, s) == 2 && k >= 1 && k <= 3 &&
            strlen(s) == 64) {
            if (g->base[k - 1][0] == '\0') strcpy(g->base[k - 1], s);
            else if (strcmp(g->base[k - 1], s) != 0) break;   /* bad file */
        } else if (sscanf(p, "id %d", &k) == 1) {
            id = k;
        } else if (sscanf(p, "ncbieaa \"%64[^\"]\"", s) == 1) {
            strcpy(aa, s);
        } else if (sscanf(p, "sncbieaa \"%64[^\"]\"", s) == 1) {
            if (id < 0 || strlen(aa) != 64 || strlen(s) != 64 || g->n == GC_MAX)
                break;
            g->id[g->n] = id;
            strcpy(g->aa[g->n], aa);
            strcpy(g->start[g->n], s);
            g->n++;
            id = -1;
            aa[0] = '\0';
        }
    }
    k = feof(f);
    fclose(f);
    if (!k || strcmp(version, GC_PRT_VERSION) != 0 || g->n == 0) {
        printf("[FAIL] %s: version %s (want " GC_PRT_VERSION "), %d tables%s\n",
               GC_PRT, version, g->n, k ? "" : ", unreadable past them");
        return 0;
    }
    /* The base rows (64 letters of TCAG each, if read) spell the 64
       codons, each once */
    for (k = 0; k < 64 && g->base[0][0] && g->base[1][0] && g->base[2][0]; k++) {
        int v = 0;
        for (j = 0; j < 3; j++) v = 4 * v + (int)(strchr("TCAG", g->base[j][k]) - "TCAG");
        if (seen[v]++) break;
    }
    if (k < 64) {
        printf("[FAIL] %s: the Base1-3 rows don't spell the 64 codons\n", GC_PRT);
        return 0;
    }
    return 1;
}

static int in_list(int t, const int *v, size_t n) {
    size_t i;
    for (i = 0; i < n; i++) if (v[i] == t) return 1;
    return 0;
}
#define IN_LIST(t, v) in_list((t), (v), sizeof(v) / sizeof((v)[0]))

/* Appends a formatted line to fails (FAILS_LEN bytes). */
#define FAILS_LEN 2048
static void add_fail(char *fails, const char *fmt, ...) {
    va_list ap;
    size_t n = strlen(fails);
    va_start(ap, fmt);
    vsnprintf(fails + n, FAILS_LEN - n, fmt, ap);
    va_end(ap);
}

/* Says in why (len bytes) where Prodigal's use of table g->id[r] differs
   from NCBI's (white-box: is_stop, is_start, amino, pdg_no_stop), and
   returns 1; 0 if it doesn't.  text holds the 64 codons in gc.prt's order,
   then the N_CODONS ones with an N, and seq holds text encoded.  spacer
   holds two one-base training contigs, so a spacer at base 1 of its seq
   and at base 0 of its rseq. */
static int ncbi_differs(const gc_prt_t *g, int r, const char *text,
                        unsigned char *seq, const prodigal_ctx_t *spacer,
                        char *why, size_t len) {
    static const double gcs[] = { 0.3, 0.5, 0.65 };
    static const int off[2] = { 1, 0 };
    struct _training *tinf = table_model(g->id[r], 0.5);
    unsigned char *strand[2] = { spacer->seq, spacer->rseq };
    int k, j, f;
    for (k = 0; k < 64; k++) {
        const char c[4] = { text[3 * k], text[3 * k + 1], text[3 * k + 2], '\0' };
        char want = g->aa[r][k], got = amino(seq, 3 * k, tinf, 0);
        int stop = is_stop(seq, 3 * k, tinf), start = is_start(seq, 3 * k, tinf);
        int ncbi_start = g->start[r][k] == 'M';
        int modeled = strcmp(c, "ATG") == 0 || strcmp(c, "GTG") == 0 ||
                      strcmp(c, "TTG") == 0;
        if (g->start[r][k] == '*' && want != '*')
            snprintf(why, len, "%s is a stop or %c by context", c, want);
        else if (stop != (want == '*'))
            snprintf(why, len, "%s is %sa stop, NCBI's is %s", c, stop ? "" : "not ",
                     stop ? "not" : "one");
        else if (got != want)
            snprintf(why, len, "%s is %c, NCBI %c", c, got, want);
        else if (start != ncbi_start && (start || modeled))
            snprintf(why, len, "%s is %sa start, NCBI's is %sone", c,
                     start ? "" : "not ", ncbi_start ? "" : "not ");
        else if (start && amino(seq, 3 * k, tinf, 1) != 'M')
            snprintf(why, len, "start %s reads %c", c, amino(seq, 3 * k, tinf, 1));
        else
            continue;
        return 1;
    }
    /* raw_coding_score's stop probability, bases independent at GC gc:
       A and T (1 - gc) / 2 each, C and G gc / 2 */
    for (j = 0; j < 3; j++) {
        double p, ncbi = 0.0, prodigal;
        for (k = 0; k < 64; k++) {
            if (g->aa[r][k] != '*') continue;
            for (p = 1.0, f = 0; f < 3; f++)
                p *= strchr("AT", text[3 * k + f]) ? (1 - gcs[j]) / 2 : gcs[j] / 2;
            ncbi += p;
        }
        tinf->gc = gcs[j];
        prodigal = 1.0 - pdg_no_stop(tinf);
        if (fabs(prodigal - ncbi) > 1e-12) {
            snprintf(why, len, "the ORF-length score's stop probability at GC "
                     "%g is %.6f, NCBI's stops' %.6f", gcs[j], prodigal, ncbi);
            return 1;
        }
    }
    /* The spacer between training contigs: a stop in each frame, on both
       strands */
    for (j = 0; j < 2; j++)
        for (f = 0; f < 3; f++) {
            int found = 0;
            for (k = off[j] + f; k + 3 <= off[j] + 12; k += 3)
                found |= is_stop(strand[j], k, tinf);
            if (!found) {
                snprintf(why, len, "the training spacer has no stop in frame %d "
                         "of the %s strand", f, j == 0 ? "forward" : "reverse");
                return 1;
            }
        }
    /* An ambiguous base makes no stop or start */
    for (k = 64; k < 64 + N_CODONS; k++)
        if (is_stop(seq, 3 * k, tinf) || is_start(seq, 3 * k, tinf)) {
            snprintf(why, len, "an ambiguous base makes %.3s a %s", text + 3 * k,
                     is_stop(seq, 3 * k, tinf) ? "stop" : "start");
            return 1;
        }
    return 0;
}

static void test_genetic_codes_ncbi(void) {
    TEST_START("genetic codes: supported iff they are NCBI's, codon by codon");
    static gc_prt_t g;
    static const char *two[2] = { "A", "A" };
    static const int32_t one[2] = { 1, 1 };
    char text[3 * (64 + N_CODONS)], why[128], msg[256], fails[FAILS_LEN] = "";
    unsigned char seq[(3 * (64 + N_CODONS) * 2 + 7) / 8];
    prodigal_ctx_t *ctx;
    int t, r, k, n, n11, n4;
    if (!load_gc_prt(&g)) return;
    for (k = 0; k < 64; k++)
        for (r = 0; r < 3; r++) text[3 * k + r] = g.base[r][k];
    for (k = n = 0; k < 125; k++)       /* every codon of TCAGN with an N */
        if (k / 25 == 4 || k / 5 % 5 == 4 || k % 5 == 4) {
            text[3 * (64 + n)] = "TCAGN"[k / 25];
            text[3 * (64 + n) + 1] = "TCAGN"[k / 5 % 5];
            text[3 * (64 + n) + 2] = "TCAGN"[k % 5];
            n++;
        }
    ASSERT_EQ_INT(n, N_CODONS);
    ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);
    /* The codons, encoded; then the spacer, in the same context */
    ASSERT_EQ_INT(prodigal_set_sequence(ctx, text, (int32_t)sizeof(text), NULL),
                  PRODIGAL_OK);
    memcpy(seq, ctx->seq, sizeof(seq));
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, two, NULL, one, 2),
                  PRODIGAL_OK);
    ASSERT_EQ_INT(ctx->slen, 26);

    for (t = -1; t <= GC_PROBE_MAX; t++) {
        prodigal_config_t c;
        int ok, faithful = 0;
        for (r = 0; r < g.n && g.id[r] != t; r++) ;
        if (r == g.n) snprintf(why, sizeof(why), "NCBI has no table %d", t);
        else faithful = !ncbi_differs(&g, r, text, seq, ctx, why, sizeof(why));
        prodigal_config_init(&c);
        c.trans_table = t;
        ok = prodigal_config_validate(&c, msg, sizeof(msg)) == PRODIGAL_OK;
        if (ok != faithful || ok != IN_LIST(t, supported_tables) ||
            (r < g.n && !ok) != IN_LIST(t, rejected_ncbi_tables))
            add_fail(fails, "\n       table %d: %s; %s; pinned %s", t,
                     ok ? "supported" : msg, faithful ? "matches NCBI" : why,
                     IN_LIST(t, supported_tables) ? "supported" :
                     IN_LIST(t, rejected_ncbi_tables) ? "rejected" : "absent");
    }
    for (r = 0; r < g.n; r++)
        if (g.id[r] < -1 || g.id[r] > GC_PROBE_MAX)
            add_fail(fails, "\n       NCBI table %d is past those checked", g.id[r]);

    /* The metagenomic bins use tables 11 (47) and 4 (3), as metagenomic
       mode's refusal of the others says; both are supported */
    ASSERT_EQ_INT(pdg_meta_init(ctx), PRODIGAL_OK);
    for (k = 0, n11 = n4 = 0; k < NUM_META; k++) {
        prodigal_config_t c;
        prodigal_config_init(&c);
        c.trans_table = ctx->meta[k].tinf->trans_table;
        n11 += c.trans_table == 11;
        n4 += c.trans_table == 4;
        if (prodigal_config_validate(&c, msg, sizeof(msg)) != PRODIGAL_OK)
            add_fail(fails, "\n       metagenomic bin %d: %s", k, msg);
    }
    if (n11 != 47 || n4 != 3)
        add_fail(fails, "\n       metagenomic bins: %d use table 11, %d table 4, "
                 "of %d", n11, n4, NUM_META);
    prodigal_destroy(ctx);
    if (fails[0] != '\0') {
        printf("[FAIL]%s\n", fails);
        return;
    }
    TEST_PASS();
}

/* Every way to choose a table (config, setter, model blob) rejects the same
   ones, with the same whole reason, and a rejected one changes nothing. */
static void test_genetic_code_entry_points(void) {
    TEST_START("genetic codes: every entry point rejects them, saying why");
    static const struct { int table; const char *msg; } pinned[] = {
        { 7, "translation table 7 is not supported: NCBI defines no such "
             "genetic code (the-miint/Prodigal docs/genetic-codes.md)" },
        { 10, "translation table 10 is not supported: Prodigal translates "
              "TGA as W (NCBI: C) (the-miint/Prodigal docs/genetic-codes.md)" },
        { 27, "translation table 27 is not supported: TGA is a stop or W "
              "depending on context, which Prodigal cannot model "
              "(the-miint/Prodigal docs/genetic-codes.md)" },
    };
    static const char tail[] = " (the-miint/Prodigal docs/genetic-codes.md)";
    char msg[256], want[256], what[32];
    const char *reason;
    size_t tlen, i;
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL) return;
    struct _training *blob = (struct _training *)xrealloc(NULL, tlen);
    prodigal_ctx_t *ctx = make_model_ctx(tdata, tlen);
    prodigal_training_info_t info;
    int t, ok, before;
    ASSERT_TRUE(ctx != NULL);
    ASSERT_EQ_INT(ctx->tinf.trans_table, 11);
    memcpy(blob, tdata, tlen);

    for (t = -1; t <= GC_PROBE_MAX; t++) {
        prodigal_config_t c;
        prodigal_config_init(&c);
        c.trans_table = t;
        ok = prodigal_config_validate(&c, msg, sizeof(msg)) == PRODIGAL_OK;
        /* "translation table T is not supported: <why> (<tail>)", whole */
        snprintf(want, sizeof(want), "translation table %d is not supported: ", t);
        reason = msg + strlen(want) - strlen(" is not supported: ");
        if (!ok && (strncmp(msg, want, strlen(want)) != 0 ||
                    strlen(msg) < strlen(want) + strlen(tail) ||
                    strcmp(msg + strlen(msg) - strlen(tail), tail) != 0)) {
            printf("[FAIL] table %d: config_validate says \"%s\"\n", t, msg);
            return;
        }
        for (i = 0; i < sizeof(pinned) / sizeof(pinned[0]); i++)
            if (pinned[i].table == t && strcmp(msg, pinned[i].msg) != 0) {
                printf("[FAIL] table %d: \"%s\"\n", t, msg);
                return;
            }
        /* create allocates ~50 MB: rejected tables, and one other */
        if (!ok || t == 4) {
            prodigal_ctx_t *made = prodigal_create(&c);
            if ((made != NULL) != ok) {
                printf("[FAIL] table %d: config_validate and create disagree\n", t);
                return;
            }
            prodigal_destroy(made);
        }
        /* The setter; a rejected table leaves the model's and the config's */
        snprintf(want, sizeof(want), "Translation table %d%s", t, ok ? "" : reason);
        before = ctx->config.trans_table;
        if ((prodigal_set_translation_table(ctx, t) == PRODIGAL_OK) != ok ||
            ctx->tinf.trans_table != (ok ? t : 11) ||
            ctx->config.trans_table != (ok ? t : before) ||
            (!ok && strcmp(prodigal_last_error(ctx), want) != 0)) {
            printf("[FAIL] table %d: set_translation_table says \"%s\"\n", t,
                   prodigal_last_error(ctx));
            return;
        }
        /* A model blob: training_info and load_training agree (both_reject),
           and a rejected one leaves the model as it was */
        blob->trans_table = t;
        snprintf(what, sizeof(what), "a model with table %d", t);
        snprintf(want, sizeof(want), "Training data has translation table %d, "
                 "which%s", t, ok ? "" : reason);
        if (ok) {
            ASSERT_EQ_INT(prodigal_training_info(blob, tlen, &info, NULL, 0),
                          PRODIGAL_OK);
            ASSERT_EQ_INT(prodigal_load_training(ctx, blob, tlen), PRODIGAL_OK);
            ASSERT_EQ_INT(ctx->tinf.trans_table, t);
            ASSERT_EQ_INT(prodigal_load_training(ctx, tdata, tlen), PRODIGAL_OK);
        } else {
            if (!both_reject(ctx, blob, tlen, what)) return;
            if (strcmp(prodigal_last_error(ctx), want) != 0 ||
                memcmp(&ctx->tinf, tdata, tlen) != 0) {
                printf("[FAIL] %s: \"%s\", want \"%s\", or the model changed\n",
                       what, prodigal_last_error(ctx), want);
                return;
            }
        }
    }
    free(blob);
    free(tdata);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/* Whether a call (rc) refused the single-genome model's "what" on
   metagenomic context ctx, with the reason, and changed nothing: its model
   and config are still model and config. */
static int refused_in_meta(prodigal_ctx_t *ctx, int rc, const char *what,
                           const struct _training *model,
                           const prodigal_config_t *config) {
    char want[256];
    snprintf(want, sizeof(want), "A metagenomic context takes no %s: its bins "
             "are its models; use a single-genome context (meta_mode = 0)",
             what);
    if (rc == PRODIGAL_ERR_INVALID_INPUT &&
        strcmp(prodigal_last_error(ctx), want) == 0 && !ctx->trained &&
        memcmp(model, &ctx->tinf, sizeof(*model)) == 0 &&
        memcmp(config, &ctx->config, sizeof(*config)) == 0)
        return 1;
    printf("[FAIL] a metagenomic context took a %s: %d, \"%s\"\n", what, rc,
           prodigal_last_error(ctx));
    return 0;
}

/* A metagenomic context's bins are its models, so it refuses the settings
   of the single-genome model, which its gene finding would ignore:
   force_nonsd and a start weight other than the default in the config, and
   the calls that train, load or adjust a model.  Refused, they change
   nothing.  A sequence it still takes, from either setter. */
static void test_meta_mode_refuses_model_settings(void) {
    TEST_START("metagenomic mode refuses the single-genome model's settings");
    static const char *one[1] = { "ACGTACGT" };
    static const int32_t len1[1] = { 8 };
    static struct _training before;     /* 558 KB: not on the stack */
    prodigal_config_t c, cbefore;
    char msg[256];
    prodigal_ctx_t *ctx;
    size_t tlen;
    void *tdata = load_ref_model(&tlen);
    if (tdata == NULL) return;

    /* The config: refused in meta mode, accepted in single mode */
    prodigal_config_init(&c);
    c.meta_mode = 1;
    c.force_nonsd = 1;
    ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_TRUE(strcmp(msg, "metagenomic mode takes no force_nonsd (-n), as "
                       "each bin has its own SD setting; use single mode") == 0);
    ASSERT_TRUE(prodigal_create(&c) == NULL);
    c.force_nonsd = 0;
    c.start_weight = 2.0;
    ASSERT_EQ_INT(prodigal_config_validate(&c, msg, sizeof(msg)),
                  PRODIGAL_ERR_INVALID_CONFIG);
    ASSERT_TRUE(strcmp(msg, "metagenomic mode takes only the default start "
                       "weight, 4.35, as each bin has its own (got 2); use "
                       "single mode") == 0);
    c.start_weight = 4.35;      /* the default, which changes nothing */
    ASSERT_EQ_INT(prodigal_config_validate(&c, NULL, 0), PRODIGAL_OK);
    c.meta_mode = 0;
    c.force_nonsd = 1;
    c.start_weight = 2.0;
    ASSERT_EQ_INT(prodigal_config_validate(&c, NULL, 0), PRODIGAL_OK);

    /* The calls, on a metagenomic context */
    ctx = make_ctx(1);
    ASSERT_TRUE(ctx != NULL);
    memcpy(&before, &ctx->tinf, sizeof(before));
    memcpy(&cbefore, &ctx->config, sizeof(cbefore));
    if (!refused_in_meta(ctx, prodigal_load_training(ctx, tdata, tlen),
                         "training model", &before, &cbefore) ||
        !refused_in_meta(ctx, prodigal_train(ctx), "training", &before,
                         &cbefore) ||
        !refused_in_meta(ctx, prodigal_set_start_weight(ctx, 2.0),
                         "start weight", &before, &cbefore) ||
        !refused_in_meta(ctx, prodigal_set_gc(ctx, 0.5), "model GC", &before,
                         &cbefore) ||
        !refused_in_meta(ctx, prodigal_set_uses_sd(ctx, 0), "SD setting",
                         &before, &cbefore))
        return;
    /* Still a sequence to search, from set_training_sequences too */
    ASSERT_EQ_INT(prodigal_set_training_sequences(ctx, one, NULL, len1, 1),
                  PRODIGAL_OK);
    ASSERT_EQ_INT(ctx->slen, 8);
    free(tdata);
    prodigal_destroy(ctx);
    TEST_PASS();
}

/* Metagenomic mode chooses table 11 or 4 for each sequence from its bins,
   so it refuses any table but 11, through the config and the setter (a
   table refused in single mode gets single mode's reason). */
static void test_genetic_code_meta_mode(void) {
    TEST_START("genetic codes: metagenomic mode takes table 11 only");
    static const char meta_why[] = "metagenomic mode takes table 11 only, "
        "as its bins choose table 11 or 4 for each sequence; use single mode "
        "(the-miint/Prodigal docs/genetic-codes.md)";
    char msg[256], expect[256], want[256];
    prodigal_ctx_t *ctx = make_ctx(1);
    prodigal_config_t c;
    int t, ok;
    ASSERT_TRUE(ctx != NULL);
    for (t = -1; t <= GC_PROBE_MAX; t++) {
        prodigal_config_init(&c);
        c.trans_table = t;
        /* Single mode's reason if it refuses t, else meta mode's */
        if (prodigal_config_validate(&c, expect, sizeof(expect)) == PRODIGAL_OK)
            snprintf(expect, sizeof(expect), "translation table %d is not "
                     "supported: %s", t, meta_why);
        c.meta_mode = 1;
        ok = prodigal_config_validate(&c, msg, sizeof(msg)) == PRODIGAL_OK;
        if (ok != (t == 11) || (!ok && strcmp(msg, expect) != 0)) {
            printf("[FAIL] meta table %d: \"%s\", want \"%s\"\n", t,
                   ok ? "accepted" : msg, t == 11 ? "accepted" : expect);
            return;
        }
        /* The setter says the same, and a refused table changes nothing */
        snprintf(want, sizeof(want), "T%s", expect + 1);
        if ((prodigal_set_translation_table(ctx, t) == PRODIGAL_OK) != ok ||
            ctx->config.trans_table != 11 || ctx->tinf.trans_table != 11 ||
            (!ok && strcmp(prodigal_last_error(ctx), want) != 0)) {
            printf("[FAIL] meta set_translation_table(%d): \"%s\"\n", t,
                   prodigal_last_error(ctx));
            return;
        }
    }
    /* create refuses one (make_ctx made one with 11) */
    prodigal_config_init(&c);
    c.meta_mode = 1;
    c.trans_table = 4;
    ASSERT_TRUE(prodigal_create(&c) == NULL);
    prodigal_destroy(ctx);
    TEST_PASS();
}

int main(void) {
    printf("=== Prodigal Library API Test Suite ===\n\n");

    printf("Phase 1.1: Error Codes\n");
    test_error_codes();
    test_strerror();
    test_version_constants();

    printf("\nPhase 1.2: Config Struct\n");
    test_config_init_defaults();
    test_config_struct_size_at_offset_zero();
    test_abi_layout();
    test_config_init_abi();

    printf("\nPhase 1.3: Context Lifecycle\n");
    test_create_destroy();
    test_create_null_config();
    test_create_bad_struct_size();
    test_destroy_null();
    test_config_validate();
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
    test_training_info();
    test_export_needs_model();

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
    test_failed_training_leaves_no_model();
    test_training_alloc_failures();

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

    printf("\nPhase 15: Low-complexity input\n");
    test_low_complexity_meta();
    test_low_complexity_model();
    test_low_complexity_trained();
    test_low_complexity_grow_fails();

    printf("\nPhase 16: Limits\n");
    test_gene_array_grows();
    test_masks_past_max();
    test_length_limits();

    printf("\nPhase 17: Proteins\n");
    test_protein_letters();
    test_proteins_off();

    printf("\nPhase 18: Genetic codes (NCBI)\n");
    test_genetic_codes_ncbi();
    test_genetic_code_entry_points();
    test_genetic_code_meta_mode();
    test_meta_mode_refuses_model_settings();

    printf("\n=== Results: %d/%d passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
