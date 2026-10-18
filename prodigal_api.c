/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Library API implementation.
*******************************************************************************/

#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _WIN32
#include <malloc.h>   /* _aligned_malloc, _aligned_free */
#endif
#include "prodigal_internal.h"

/*******************************************************************************
    Validation helpers
*******************************************************************************/

static int is_valid_trans_table(int tt) {
    if (tt < 1 || tt > 25) return 0;
    if (tt == 7 || tt == 8) return 0;
    if (tt >= 17 && tt <= 20) return 0;
    return 1;
}

/*******************************************************************************
    Error reporting
*******************************************************************************/

const char *prodigal_version_string(void) {
    return PRODIGAL_VERSION_STRING;
}

const char *prodigal_strerror(int error_code) {
    switch (error_code) {
        case PRODIGAL_OK:                  return "Success";
        case PRODIGAL_ERR_NOMEM:           return "Out of memory";
        case PRODIGAL_ERR_INVALID_CONFIG:  return "Invalid configuration";
        case PRODIGAL_ERR_INVALID_INPUT:   return "Invalid input";
        case PRODIGAL_ERR_INTERNAL:        return "Internal error";
        case PRODIGAL_ERR_SEQ_TOO_SHORT:   return "Sequence too short";
        case PRODIGAL_ERR_CANCELLED:       return "Cancelled";
        default:                           return "Unknown error";
    }
}

const char *prodigal_last_error(const prodigal_ctx_t *ctx) {
    if (ctx == NULL) return "NULL context";
    return ctx->error_msg;
}

/*******************************************************************************
    Config initialization
*******************************************************************************/

void prodigal_config_init(prodigal_config_t *config) {
    memset(config, 0, sizeof(*config));
    config->struct_size = sizeof(prodigal_config_t);
    config->trans_table = 11;
    config->start_weight = 4.35;
}

/*******************************************************************************
    Context lifecycle
*******************************************************************************/

static void free_meta(prodigal_ctx_t *ctx) {
    int i;
    if (ctx->meta == NULL) return;
    for (i = 0; i < NUM_META; i++) pdg_free(ctx, ctx->meta[i].tinf);
    pdg_free(ctx, ctx->meta);
    ctx->meta = NULL;
}

/* A fresh model with the context's training parameters, which the config
   holds (the setters update it). */
static void reset_model(prodigal_ctx_t *ctx) {
    memset(&ctx->tinf, 0, sizeof(ctx->tinf));
    ctx->tinf.trans_table = ctx->config.trans_table;
    ctx->tinf.st_wt = ctx->config.start_weight;
}

prodigal_ctx_t *prodigal_create(const prodigal_config_t *config) {
    if (config == NULL) return NULL;
    if (config->struct_size != sizeof(prodigal_config_t)) return NULL;
    if (!is_valid_trans_table(config->trans_table)) return NULL;

    prodigal_ctx_t *ctx = (prodigal_ctx_t *)malloc(sizeof(prodigal_ctx_t));
    if (ctx == NULL) return NULL;
    memset(ctx, 0, sizeof(*ctx));

    /* Snapshot config */
    ctx->config = *config;
    reset_model(ctx);

    /* Allocate sequence buffers */
    ctx->seq = (unsigned char *)pdg_alloc(ctx, MAX_SEQ / 4 * sizeof(unsigned char));
    ctx->rseq = (unsigned char *)pdg_alloc(ctx, MAX_SEQ / 4 * sizeof(unsigned char));
    ctx->useq = (unsigned char *)pdg_alloc(ctx, MAX_SEQ / 8 * sizeof(unsigned char));
    if (ctx->seq == NULL || ctx->rseq == NULL || ctx->useq == NULL) {
        prodigal_destroy(ctx);
        return NULL;
    }
    memset(ctx->seq, 0, MAX_SEQ / 4 * sizeof(unsigned char));
    memset(ctx->rseq, 0, MAX_SEQ / 4 * sizeof(unsigned char));
    memset(ctx->useq, 0, MAX_SEQ / 8 * sizeof(unsigned char));

    /* Allocate node and gene arrays */
    ctx->nodes = (struct _node *)pdg_alloc(ctx, STT_NOD * sizeof(struct _node));
    ctx->genes = (struct _gene *)pdg_alloc(ctx, MAX_GENES * sizeof(struct _gene));
    if (ctx->nodes == NULL || ctx->genes == NULL) {
        prodigal_destroy(ctx);
        return NULL;
    }
    memset(ctx->nodes, 0, STT_NOD * sizeof(struct _node));
    memset(ctx->genes, 0, MAX_GENES * sizeof(struct _gene));

    return ctx;
}

void prodigal_destroy(prodigal_ctx_t *ctx) {
    if (ctx == NULL) return;

    pdg_free(ctx, ctx->seq);
    pdg_free(ctx, ctx->rseq);
    pdg_free(ctx, ctx->useq);
    pdg_free(ctx, ctx->nodes);
    pdg_free(ctx, ctx->genes);

    free_meta(ctx);

    free(ctx);  /* Context itself always uses system malloc */
}

/*******************************************************************************
    Sequence encoding helpers
*******************************************************************************/

/* Encode a single base into the bitmap at position *bctr, advance counters.
   Returns 1 if base is G or C (for GC counting), 0 otherwise. */
static int encode_base(unsigned char *seq, unsigned char *useq,
                       int *bctr, int *len, char ch) {
    int gc = 0;
    if (ch == 'g' || ch == 'G') {
        set(seq, *bctr);
        gc = 1;
    }
    else if (ch == 't' || ch == 'T') {
        set(seq, *bctr);
        set(seq, *bctr + 1);
    }
    else if (ch == 'c' || ch == 'C') {
        set(seq, *bctr + 1);
        gc = 1;
    }
    else if (ch != 'a' && ch != 'A') {
        /* Ambiguous base: encode as C, mark in useq */
        set(seq, *bctr + 1);
        set(useq, *len);
    }
    /* A = nothing set (00) */
    *bctr += 2;
    (*len)++;
    return gc;
}

/* Insert TTAATTAATTAA (12 bases) as stop codons in all 6 frames */
static void encode_stop_spacer(unsigned char *seq, int *bctr, int *len) {
    int i;
    for (i = 0; i < 12; i++) {
        if (i % 4 == 0 || i % 4 == 1) {
            set(seq, *bctr);
            set(seq, *bctr + 1);
        }
        *bctr += 2;
        (*len)++;
    }
}

void pdg_clear_sequence(prodigal_ctx_t *ctx) {
    memset(ctx->seq, 0, (ctx->slen / 4 + 1) * sizeof(unsigned char));
    memset(ctx->rseq, 0, (ctx->slen / 4 + 1) * sizeof(unsigned char));
    memset(ctx->useq, 0, (ctx->slen / 8 + 1) * sizeof(unsigned char));
    ctx->slen = 0;
    ctx->nmask = 0;
    ctx->gc = 0.0;
    ctx->cur_header[0] = '\0';
    ctx->error_msg[0] = '\0';
    ctx->error_code = PRODIGAL_OK;
}

/*******************************************************************************
    Sequence input
*******************************************************************************/

int prodigal_set_sequence(prodigal_ctx_t *ctx, const char *seq, int32_t len,
                          const char *header) {
    int bctr = 0, slen = 0, gc_cont = 0;
    int mask_beg = -1;
    int32_t i;

    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (seq == NULL || len <= 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Sequence is NULL or length <= 0");
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    if (len >= MAX_SEQ) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Sequence length %d exceeds maximum %d", len, MAX_SEQ);
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    pdg_clear_sequence(ctx);

    /* Encode each base */
    for (i = 0; i < len; i++) {
        char ch = seq[i];
        if (ch < 'A' || ch > 'z') continue;

        /* Masking logic */
        if (ctx->config.mask_regions) {
            if (mask_beg != -1 && ch != 'N' && ch != 'n') {
                if (slen - mask_beg >= MASK_SIZE) {
                    if (ctx->nmask < MAX_MASKS) {
                        ctx->mlist[ctx->nmask].begin = mask_beg;
                        ctx->mlist[ctx->nmask].end = slen - 1;
                        ctx->nmask++;
                    }
                }
                mask_beg = -1;
            }
            if (mask_beg == -1 && (ch == 'N' || ch == 'n'))
                mask_beg = slen;
        }

        gc_cont += encode_base(ctx->seq, ctx->useq, &bctr, &slen, ch);
    }

    if (slen == 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "No valid bases found in sequence");
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    ctx->slen = slen;
    ctx->gc = (double)gc_cont / (double)slen;

    /* Compute reverse complement */
    rcom_seq(ctx->seq, ctx->rseq, ctx->useq, slen);

    /* Store header */
    if (header != NULL) {
        strncpy(ctx->cur_header, header, MAX_LINE - 1);
        ctx->cur_header[MAX_LINE - 1] = '\0';
    }

    return PRODIGAL_OK;
}

int prodigal_set_training_sequences(prodigal_ctx_t *ctx, const char **seqs,
                                    const char **headers, const int32_t *lens,
                                    int32_t n_seqs) {
    int bctr = 0, slen = 0, gc_cont = 0;
    int mask_beg = -1;
    int32_t s, i;

    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (seqs == NULL || lens == NULL || n_seqs <= 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Invalid training sequence arguments");
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    pdg_clear_sequence(ctx);

    for (s = 0; s < n_seqs; s++) {
        /* Insert stop spacer between sequences (and after last) */
        if (s > 0) {
            encode_stop_spacer(ctx->seq, &bctr, &slen);
        }

        if (seqs[s] == NULL || lens[s] <= 0) continue;

        for (i = 0; i < lens[s]; i++) {
            char ch = seqs[s][i];
            if (ch < 'A' || ch > 'z') continue;

            /* Masking logic */
            if (ctx->config.mask_regions) {
                if (mask_beg != -1 && ch != 'N' && ch != 'n') {
                    if (slen - mask_beg >= MASK_SIZE) {
                        if (ctx->nmask < MAX_MASKS) {
                            ctx->mlist[ctx->nmask].begin = mask_beg;
                            ctx->mlist[ctx->nmask].end = slen - 1;
                            ctx->nmask++;
                        }
                    }
                    mask_beg = -1;
                }
                if (mask_beg == -1 && (ch == 'N' || ch == 'n'))
                    mask_beg = slen;
            }

            gc_cont += encode_base(ctx->seq, ctx->useq, &bctr, &slen, ch);

            if (slen + MAX_LINE >= MAX_SEQ) break;
        }
        if (slen + MAX_LINE >= MAX_SEQ) break;
    }

    /* Trailing stop spacer if multiple sequences */
    if (n_seqs > 1) {
        encode_stop_spacer(ctx->seq, &bctr, &slen);
    }

    if (slen == 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "No valid bases found in training sequences");
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    ctx->slen = slen;
    ctx->gc = (double)gc_cont / (double)slen;
    rcom_seq(ctx->seq, ctx->rseq, ctx->useq, slen);

    /* Store first header */
    if (headers != NULL && headers[0] != NULL) {
        strncpy(ctx->cur_header, headers[0], MAX_LINE - 1);
        ctx->cur_header[MAX_LINE - 1] = '\0';
    }

    return PRODIGAL_OK;
}

int prodigal_get_seq_info(const prodigal_ctx_t *ctx, prodigal_seq_info_t *info) {
    if (ctx == NULL || info == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    info->length = ctx->slen;
    info->gc_content = ctx->gc;
    return PRODIGAL_OK;
}

/*******************************************************************************
    Nodes
*******************************************************************************/

/* Make room for the nodes of the current sequence. A new array is zeroed. */
static int reserve_nodes(prodigal_ctx_t *ctx) {
    struct _node *nodes;
    if (ctx->slen <= ctx->max_slen || ctx->slen <= STT_NOD * 8) return PRODIGAL_OK;
    nodes = (struct _node *)pdg_calloc(ctx,
        (size_t)(ctx->slen / 8) * sizeof(struct _node));
    if (nodes == NULL) {
        /* Upstream's message, which the CLI prints */
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM, "Realloc failed on nodes");
        return PRODIGAL_ERR_NOMEM;
    }
    pdg_free(ctx, ctx->nodes);
    ctx->nodes = nodes;
    ctx->nn = 0;
    ctx->max_slen = ctx->slen;
    return PRODIGAL_OK;
}

/* The sorted nodes of the current sequence under model tinf, in ctx->nodes.
   add_nodes and the steps after it expect zeroed nodes (add_nodes sets edge
   only when it is 1, and a stale elim drops a gene), so the previous list,
   ctx->nodes[0, ctx->nn), is cleared first; everything past it is zero. */
static int build_nodes(prodigal_ctx_t *ctx, struct _training *tinf) {
    int rc = reserve_nodes(ctx);
    if (rc != PRODIGAL_OK) return rc;
    memset(ctx->nodes, 0, ctx->nn * sizeof(struct _node));
    ctx->nn = add_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                        ctx->config.closed_ends, ctx->mlist, ctx->nmask, tinf);
    qsort(ctx->nodes, ctx->nn, sizeof(struct _node), &compare_nodes);
    return PRODIGAL_OK;
}

/*******************************************************************************
    Training
*******************************************************************************/

int pdg_train(prodigal_ctx_t *ctx, pdg_train_report_fn report) {
    int *gc_frame;
    int ipath, rc;

    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (ctx->slen == 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "No sequence loaded for training");
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    if (ctx->slen < 20000) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_SEQ_TOO_SHORT,
                      "Sequence must be >= 20000 bp for training (got %d)",
                      ctx->slen);
        return PRODIGAL_ERR_SEQ_TOO_SHORT;
    }
    /* Start from a fresh model, as the CLI does, with the training
       sequence's GC (the steps below read it). Until training succeeds
       there is no model. */
    reset_model(ctx);
    ctx->tinf.gc = ctx->gc;
    ctx->trained = 0;

#define REPORT(stage, done) do { \
    if (report != NULL) report(ctx, (stage), (done)); \
} while (0)

    REPORT(PDG_TRAIN_NODES, 0);
    rc = build_nodes(ctx, &ctx->tinf);
    if (rc != PRODIGAL_OK) return rc;
    REPORT(PDG_TRAIN_NODES, 1);

    /* GC frame bias */
    REPORT(PDG_TRAIN_GC_BIAS, 0);
    gc_frame = calc_most_gc_frame(ctx->seq, ctx->slen);
    if (gc_frame == NULL) {
        /* Upstream's message, which the CLI prints */
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM, "Malloc failed on gc frame plot");
        return PRODIGAL_ERR_NOMEM;
    }
    record_gc_bias(gc_frame, ctx->nodes, ctx->nn, &ctx->tinf);
    free(gc_frame);
    REPORT(PDG_TRAIN_GC_BIAS, 1);

    /* Initial DP with GC bias only */
    REPORT(PDG_TRAIN_INITIAL_GENES, 0);
    record_overlapping_starts(ctx->nodes, ctx->nn, &ctx->tinf, 0);
    ipath = dprog(ctx->nodes, ctx->nn, &ctx->tinf, 0);
    REPORT(PDG_TRAIN_INITIAL_GENES, 1);

    /* Dicodon statistics and coding scores */
    REPORT(PDG_TRAIN_CODING, 0);
    calc_dicodon_gene(&ctx->tinf, ctx->seq, ctx->rseq, ctx->slen,
                      ctx->nodes, ipath);
    raw_coding_score(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                     ctx->nn, &ctx->tinf);
    REPORT(PDG_TRAIN_CODING, 1);

    /* RBS and start training */
    REPORT(PDG_TRAIN_STARTS, 0);
    rbs_score(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, ctx->nn, &ctx->tinf);
    train_starts_sd(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, ctx->nn,
                    &ctx->tinf);
    determine_sd_usage(&ctx->tinf);
    if (ctx->config.force_nonsd) ctx->tinf.uses_sd = 0;
    if (ctx->tinf.uses_sd == 0)
        train_starts_nonsd(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                           ctx->nn, &ctx->tinf);
    ctx->trained = 1;
    REPORT(PDG_TRAIN_STARTS, 1);
#undef REPORT

    return PRODIGAL_OK;
}

/* prodigal_train's report: the log callback's messages. */
static void log_training(const prodigal_ctx_t *ctx, int stage, int done) {
    static const char *const what[] = {
        "Finding all potential starts and stops...",
        "Looking for GC bias in different frames...",
        "Building initial gene set...",
        "Creating coding model and scoring nodes...",
        "Examining upstream regions and training starts..."
    };
    if (!done) PDG_LOG(ctx, "%s", what[stage]);
    else if (stage == PDG_TRAIN_NODES) PDG_LOG(ctx, "%d nodes found", ctx->nn);
    else if (stage == PDG_TRAIN_STARTS)
        PDG_LOG(ctx, "Training complete (GC=%.2f, uses_sd=%d)",
                ctx->tinf.gc, ctx->tinf.uses_sd);
}

int prodigal_train(prodigal_ctx_t *ctx) {
    return pdg_train(ctx, log_training);
}

int prodigal_load_training(prodigal_ctx_t *ctx, const void *data, size_t len) {
    int table;
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (data == NULL || len != sizeof(struct _training)) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Invalid training data (expected %zu bytes, got %zu)",
                      sizeof(struct _training), len);
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    /* Check before copying, so a rejected blob changes nothing (the blob
       need not be aligned) */
    memcpy(&table, (const char *)data + offsetof(struct _training, trans_table),
           sizeof(table));
    if (!is_valid_trans_table(table)) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Training data has invalid translation table %d", table);
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    memcpy(&ctx->tinf, data, sizeof(struct _training));
    ctx->trained = 1;
    return PRODIGAL_OK;
}

int prodigal_export_training(const prodigal_ctx_t *ctx, void **data_out,
                             size_t *len_out) {
    void *buf;
    if (ctx == NULL || data_out == NULL || len_out == NULL)
        return PRODIGAL_ERR_INVALID_INPUT;

    buf = malloc(sizeof(struct _training));
    if (buf == NULL) return PRODIGAL_ERR_NOMEM;

    memcpy(buf, &ctx->tinf, sizeof(struct _training));
    *data_out = buf;
    *len_out = sizeof(struct _training);
    return PRODIGAL_OK;
}

/*******************************************************************************
    Training parameter setters
*******************************************************************************/

int prodigal_set_translation_table(prodigal_ctx_t *ctx, int table) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (!is_valid_trans_table(table)) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Invalid translation table %d", table);
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    ctx->config.trans_table = ctx->tinf.trans_table = table;
    return PRODIGAL_OK;
}

int prodigal_set_start_weight(prodigal_ctx_t *ctx, double weight) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (weight <= 0.0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Start weight must be > 0 (got %.2f)", weight);
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    ctx->config.start_weight = ctx->tinf.st_wt = weight;
    return PRODIGAL_OK;
}

int prodigal_set_gc(prodigal_ctx_t *ctx, double gc) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (gc < 0.0 || gc > 1.0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "GC must be in [0, 1] (got %.4f)", gc);
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    ctx->tinf.gc = gc;
    return PRODIGAL_OK;
}

int prodigal_set_uses_sd(prodigal_ctx_t *ctx, int uses_sd) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    ctx->tinf.uses_sd = (uses_sd != 0) ? 1 : 0;
    return PRODIGAL_OK;
}

/*******************************************************************************
    The pipeline shared with the CLI (prodigal_internal.h)
*******************************************************************************/

int pdg_score_nodes(prodigal_ctx_t *ctx, struct _training *tinf) {
    int rc = build_nodes(ctx, tinf);
    if (rc == PRODIGAL_OK)
        score_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, ctx->nn,
                    tinf, ctx->config.closed_ends, ctx->config.meta_mode != 0);
    return rc;
}

/* The genes on the DP path that ends at node ipath (-1: none), into
   ctx->genes with their attributes under model tinf; returns how many. */
static int keep_genes(prodigal_ctx_t *ctx, int ipath, struct _training *tinf,
                      int seqnum) {
    int ng;
    eliminate_bad_genes(ctx->nodes, ipath, tinf);
    ng = add_genes(ctx->genes, ctx->nodes, ipath);
    tweak_final_starts(ctx->genes, ng, ctx->nodes, ctx->nn, tinf);
    record_gene_data(ctx->genes, ng, ctx->nodes, tinf, seqnum);
    return ng;
}

int pdg_single_genes(prodigal_ctx_t *ctx, int seqnum) {
    int ipath;
    record_overlapping_starts(ctx->nodes, ctx->nn, &ctx->tinf, 1);
    ipath = dprog(ctx->nodes, ctx->nn, &ctx->tinf, 1);
    return keep_genes(ctx, ipath, &ctx->tinf, seqnum);
}

int pdg_meta_init(prodigal_ctx_t *ctx) {
    int i;
    if (ctx->meta != NULL) return PRODIGAL_OK;
    ctx->meta = (struct _metagenomic_bin *)pdg_calloc(ctx,
        NUM_META * sizeof(struct _metagenomic_bin));
    if (ctx->meta == NULL) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                      "Malloc failed on metagenomic bins.");
        return PRODIGAL_ERR_NOMEM;
    }
    /* initialize_metagenomic_bins writes every byte of each model */
    for (i = 0; i < NUM_META; i++) {
        strcpy(ctx->meta[i].desc, "None");
        ctx->meta[i].tinf = (struct _training *)pdg_alloc(ctx,
            sizeof(struct _training));
        if (ctx->meta[i].tinf == NULL) {
            free_meta(ctx);
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                          "Malloc failed on training structure.");
            return PRODIGAL_ERR_NOMEM;
        }
    }
    initialize_metagenomic_bins(ctx->meta);
    return PRODIGAL_OK;
}

int pdg_meta(prodigal_ctx_t *ctx, int seqnum, pdg_meta_result_t *res) {
    double low, high, max_score = -100.0, path_score;
    struct _training *tinf;
    int i, ipath, rc;

    res->n_genes = 0;
    res->bin = 0;
    res->n_nodes = 0;
    rc = pdg_meta_init(ctx);
    if (rc != PRODIGAL_OK) return rc;

    /* The bins whose GC is near the sequence's; every GC has some */
    low = 0.88495 * ctx->gc - 0.0102337;
    if (low > 0.65) low = 0.65;
    high = 0.86596 * ctx->gc + 0.1131991;
    if (high < 0.35) high = 0.35;

    for (i = 0; i < NUM_META; i++) {
        tinf = ctx->meta[i].tinf;
        if (ctx->config.progress_callback &&
            ctx->config.progress_callback("metagenomic scoring",
                (double)i / NUM_META, ctx->config.progress_user_data)) {
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_CANCELLED, "Cancelled");
            return PRODIGAL_ERR_CANCELLED;
        }

        /* The nodes are rebuilt only when the translation table changes,
           and before the GC window check, so a bin scores the nodes as the
           bins before it left them (docs/upstream-behavior.md, quirk 1). */
        if (i == 0 || tinf->trans_table != ctx->meta[i - 1].tinf->trans_table) {
            rc = build_nodes(ctx, tinf);
            if (rc != PRODIGAL_OK) return rc;
        }
        if (tinf->gc < low || tinf->gc > high) continue;

        reset_node_scores(ctx->nodes, ctx->nn);
        score_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes, ctx->nn,
                    tinf, ctx->config.closed_ends, 1);
        record_overlapping_starts(ctx->nodes, ctx->nn, tinf, 1);
        ipath = dprog(ctx->nodes, ctx->nn, tinf, 1);

        /* A bin without a gene path (ipath -1) scores 0.0, the empty path
           (docs/upstream-behavior.md). */
        path_score = (ipath >= 0) ? ctx->nodes[ipath].score : 0.0;
        if (path_score > max_score) {
            max_score = path_score;
            res->bin = i;
            res->n_nodes = ctx->nn;
            res->n_genes = keep_genes(ctx, ipath, tinf, seqnum);
        }
    }
    return PRODIGAL_OK;
}

/*******************************************************************************
    SOA allocation and extraction
*******************************************************************************/

#define ALIGN16(x) (((x) + 15) & ~(size_t)15)

/* The output for genes[0, ng), as record_gene_data captured them. */
static prodigal_genes_soa_t *extract_soa(const struct _gene *genes, int ng) {
    prodigal_genes_soa_t *soa;
    size_t n = (size_t)(ng > 0 ? ng : 1);  /* at least 1 to avoid zero alloc */
    size_t offset = 0;
    char *base;
    struct _rbs_text *text;
    int i;

    /* Calculate total size with 16-byte alignment */
    size_t sz_begin    = ALIGN16(n * sizeof(int32_t));
    size_t sz_end      = ALIGN16(n * sizeof(int32_t));
    size_t sz_strand   = ALIGN16(n * sizeof(int32_t));
    size_t sz_pleft    = ALIGN16(n * sizeof(int32_t));
    size_t sz_pright   = ALIGN16(n * sizeof(int32_t));
    size_t sz_stype    = ALIGN16(n * sizeof(int32_t));
    size_t sz_cscore   = ALIGN16(n * sizeof(double));
    size_t sz_sscore   = ALIGN16(n * sizeof(double));
    size_t sz_rscore   = ALIGN16(n * sizeof(double));
    size_t sz_uscore   = ALIGN16(n * sizeof(double));
    size_t sz_tscore   = ALIGN16(n * sizeof(double));
    size_t sz_conf     = ALIGN16(n * sizeof(double));
    size_t sz_gc       = ALIGN16(n * sizeof(double));
    size_t sz_motif    = ALIGN16(n * sizeof(const char *));
    size_t sz_spacer   = ALIGN16(n * sizeof(const char *));
    size_t sz_text     = ALIGN16(n * sizeof(struct _rbs_text));

    size_t total = sz_begin + sz_end + sz_strand + sz_pleft + sz_pright +
                   sz_stype + sz_cscore + sz_sscore + sz_rscore + sz_uscore +
                   sz_tscore + sz_conf + sz_gc + sz_motif + sz_spacer +
                   sz_text;

    soa = (prodigal_genes_soa_t *)malloc(sizeof(prodigal_genes_soa_t));
    if (soa == NULL) return NULL;
    memset(soa, 0, sizeof(*soa));
    soa->n_genes = ng;

    if (ng == 0) {
        soa->_base = NULL;
        return soa;
    }

    /* Single aligned allocation — portable across POSIX, Windows, WASM */
#if defined(_WIN32)
    base = (char *)_aligned_malloc(total, 16);
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__APPLE__)
    base = (char *)aligned_alloc(16, total);
#else
    {
        int pma_rv = posix_memalign((void **)&base, 16, total);
        if (pma_rv != 0) base = NULL;
    }
#endif
    if (base == NULL) {
        free(soa);
        return NULL;
    }
    soa->_base = base;   /* every element below ng is written */

    /* Carve out sub-arrays */
    offset = 0;
    soa->begin        = (int32_t *)(base + offset); offset += sz_begin;
    soa->end          = (int32_t *)(base + offset); offset += sz_end;
    soa->strand       = (int32_t *)(base + offset); offset += sz_strand;
    soa->partial_left = (int32_t *)(base + offset); offset += sz_pleft;
    soa->partial_right= (int32_t *)(base + offset); offset += sz_pright;
    soa->start_type   = (int32_t *)(base + offset); offset += sz_stype;
    soa->cscore       = (double *)(base + offset);  offset += sz_cscore;
    soa->sscore       = (double *)(base + offset);  offset += sz_sscore;
    soa->rscore       = (double *)(base + offset);  offset += sz_rscore;
    soa->uscore       = (double *)(base + offset);  offset += sz_uscore;
    soa->tscore       = (double *)(base + offset);  offset += sz_tscore;
    soa->confidence   = (double *)(base + offset);  offset += sz_conf;
    soa->gc_cont      = (double *)(base + offset);  offset += sz_gc;
    soa->rbs_motif    = (const char **)(base + offset); offset += sz_motif;
    soa->rbs_spacer   = (const char **)(base + offset); offset += sz_spacer;
    text              = (struct _rbs_text *)(base + offset);

    for (i = 0; i < ng; i++) {
        const struct _gene *g = &genes[i];
        soa->begin[i] = g->begin;
        soa->end[i] = g->end;
        soa->strand[i] = g->strand;
        soa->partial_left[i] = g->partial_left;
        soa->partial_right[i] = g->partial_right;
        soa->start_type[i] = g->start_type;
        soa->cscore[i] = g->cscore;
        soa->sscore[i] = g->sscore;
        soa->rscore[i] = g->rscore;
        soa->uscore[i] = g->uscore;
        soa->tscore[i] = g->tscore;
        soa->confidence[i] = g->confidence;
        soa->gc_cont[i] = g->gc_cont;
        text[i] = g->rbs;
        soa->rbs_motif[i] = text[i].motif;
        soa->rbs_spacer[i] = text[i].spacer;
    }

    return soa;
}

/*******************************************************************************
    Gene finding — public API
*******************************************************************************/

/* bin is the metagenomic bin that called the genes, -1 in single mode. */
static void fill_stats(const prodigal_ctx_t *ctx, int ng, int nn, int bin,
                       prodigal_stats_t *stats) {
    const struct _training *tinf = (bin >= 0) ? ctx->meta[bin].tinf : &ctx->tinf;
    if (stats == NULL) return;
    memset(stats, 0, sizeof(*stats));
    stats->n_genes = ng;
    stats->n_nodes = nn;
    stats->gc_content = ctx->gc;
    stats->translation_table = tinf->trans_table;
    stats->uses_sd = tinf->uses_sd;
    stats->best_meta_bin = bin;
    if (bin >= 0)
        snprintf(stats->best_meta_desc, sizeof(stats->best_meta_desc),
                 "%s", ctx->meta[bin].desc);
}

int prodigal_find_genes(prodigal_ctx_t *ctx, prodigal_genes_soa_t **genes_out,
                        prodigal_stats_t *stats_out) {
    int rc, ng, nn, bin = -1;

    if (ctx == NULL || genes_out == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    *genes_out = NULL;
    if (!ctx->config.meta_mode && !ctx->trained) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "No training data loaded (call prodigal_train or "
                      "prodigal_load_training first)");
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    if (ctx->slen == 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT, "No sequence loaded");
        return PRODIGAL_ERR_INVALID_INPUT;
    }

    if (ctx->config.meta_mode) {
        pdg_meta_result_t res;
        rc = pdg_meta(ctx, 1, &res);
        if (rc != PRODIGAL_OK) return rc;
        ng = res.n_genes;
        nn = res.n_nodes;
        bin = res.bin;
    } else {
        rc = pdg_score_nodes(ctx, &ctx->tinf);
        if (rc != PRODIGAL_OK) return rc;
        ng = pdg_single_genes(ctx, 1);
        nn = ctx->nn;
    }

    *genes_out = extract_soa(ctx->genes, ng);
    if (*genes_out == NULL) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                      "Failed to allocate SOA output");
        return PRODIGAL_ERR_NOMEM;
    }
    fill_stats(ctx, ng, nn, bin, stats_out);
    return PRODIGAL_OK;
}

int prodigal_find_genes_aos(prodigal_ctx_t *ctx, prodigal_genes_t **genes_out,
                            prodigal_stats_t *stats_out) {
    prodigal_genes_soa_t *soa = NULL;
    prodigal_genes_t *aos;
    int rc, i;

    if (genes_out == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    *genes_out = NULL;
    rc = prodigal_find_genes(ctx, &soa, stats_out);
    if (rc != PRODIGAL_OK) return rc;

    /* Convert SOA to AOS */
    aos = (prodigal_genes_t *)malloc(sizeof(prodigal_genes_t));
    if (aos == NULL) {
        prodigal_genes_free(soa);
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM, "Failed to allocate AOS output");
        return PRODIGAL_ERR_NOMEM;
    }
    memset(aos, 0, sizeof(*aos));
    aos->n_genes = soa->n_genes;

    if (soa->n_genes > 0) {
        struct _rbs_text *text;
        /* The genes, then their motif and spacer text */
        aos->_base = malloc(soa->n_genes * (sizeof(prodigal_gene_t) +
                                            sizeof(struct _rbs_text)));
        if (aos->_base == NULL) {
            free(aos);
            prodigal_genes_free(soa);
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                          "Failed to allocate AOS output");
            return PRODIGAL_ERR_NOMEM;
        }
        aos->genes = (prodigal_gene_t *)aos->_base;
        text = (struct _rbs_text *)(aos->genes + soa->n_genes);

        for (i = 0; i < soa->n_genes; i++) {
            aos->genes[i].begin = soa->begin[i];
            aos->genes[i].end = soa->end[i];
            aos->genes[i].strand = soa->strand[i];
            aos->genes[i].partial_left = soa->partial_left[i];
            aos->genes[i].partial_right = soa->partial_right[i];
            aos->genes[i].start_type = soa->start_type[i];
            aos->genes[i].cscore = soa->cscore[i];
            aos->genes[i].sscore = soa->sscore[i];
            aos->genes[i].rscore = soa->rscore[i];
            aos->genes[i].uscore = soa->uscore[i];
            aos->genes[i].tscore = soa->tscore[i];
            aos->genes[i].confidence = soa->confidence[i];
            aos->genes[i].gc_cont = soa->gc_cont[i];
            /* The SOA's strings came from a struct _rbs_text, so they fit */
            aos->genes[i].rbs_motif = strcpy(text[i].motif, soa->rbs_motif[i]);
            aos->genes[i].rbs_spacer = strcpy(text[i].spacer, soa->rbs_spacer[i]);
        }
    }

    prodigal_genes_free(soa);
    *genes_out = aos;
    return PRODIGAL_OK;
}

/*******************************************************************************
    Output cleanup
*******************************************************************************/

static void free_aligned(void *ptr) {
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

void prodigal_genes_free(prodigal_genes_soa_t *genes) {
    if (genes == NULL) return;
    free_aligned(genes->_base);
    free(genes);
}

void prodigal_genes_aos_free(prodigal_genes_t *genes) {
    if (genes == NULL) return;
    free(genes->_base);
    free(genes);
}
