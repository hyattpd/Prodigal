/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Library API implementation.
*******************************************************************************/

#include <stddef.h>
#include <stdint.h>
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

static int is_valid_start_weight(double w) { return w > 0.0 && isfinite(w); }
static int is_valid_gc(double gc) { return gc >= 0.0 && gc <= 1.0; } /* not NaN */

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

int prodigal_config_init_v(prodigal_config_t *config, int abi_version,
                           size_t struct_size) {
    /* Nothing is written for a caller with another ABI */
    if (config == NULL || abi_version != PRODIGAL_ABI_VERSION ||
        struct_size != sizeof(prodigal_config_t))
        return PRODIGAL_ERR_INVALID_CONFIG;
    memset(config, 0, sizeof(*config));
    config->struct_size = sizeof(prodigal_config_t);
    config->trans_table = 11;
    config->start_weight = 4.35;
    return PRODIGAL_OK;
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

int prodigal_config_validate(const prodigal_config_t *config, char *msg,
                             size_t msg_len) {
    char buf[256] = "";
    const char *why = NULL;
    if (config == NULL)
        snprintf(buf, sizeof(buf), "config is NULL");
    else if (config->struct_size != sizeof(prodigal_config_t))
        snprintf(buf, sizeof(buf), "struct_size is %zu, this library's is %zu "
                 "(header and library differ)", config->struct_size,
                 sizeof(prodigal_config_t));
    else if ((why = pdg_trans_table_unsupported(config->trans_table)) != NULL)
        snprintf(buf, sizeof(buf), "translation table %d"
                 PDG_TABLE_UNSUPPORTED, config->trans_table, why);
    else if (!is_valid_start_weight(config->start_weight))
        snprintf(buf, sizeof(buf), "start weight must be finite and > 0 "
                 "(got %g)", config->start_weight);
    else if (config->free_fn != NULL && config->alloc_fn == NULL)
        snprintf(buf, sizeof(buf), "free_fn without alloc_fn");
    if (msg != NULL && msg_len > 0) snprintf(msg, msg_len, "%s", buf);
    return buf[0] == '\0' ? PRODIGAL_OK : PRODIGAL_ERR_INVALID_CONFIG;
}

prodigal_ctx_t *prodigal_create(const prodigal_config_t *config) {
    if (prodigal_config_validate(config, NULL, 0) != PRODIGAL_OK) return NULL;

    prodigal_ctx_t *ctx = (prodigal_ctx_t *)malloc(sizeof(prodigal_ctx_t));
    if (ctx == NULL) return NULL;
    memset(ctx, 0, sizeof(*ctx));

    /* Snapshot config */
    ctx->config = *config;
    reset_model(ctx);

    /* Allocate sequence buffers */
    ctx->seq = (unsigned char *)pdg_alloc(ctx, PDG_SEQ_CAP / 4 * sizeof(unsigned char));
    ctx->rseq = (unsigned char *)pdg_alloc(ctx, PDG_SEQ_CAP / 4 * sizeof(unsigned char));
    ctx->useq = (unsigned char *)pdg_alloc(ctx, PDG_SEQ_CAP / 8 * sizeof(unsigned char));
    if (ctx->seq == NULL || ctx->rseq == NULL || ctx->useq == NULL) {
        prodigal_destroy(ctx);
        return NULL;
    }
    memset(ctx->seq, 0, PDG_SEQ_CAP / 4 * sizeof(unsigned char));
    memset(ctx->rseq, 0, PDG_SEQ_CAP / 4 * sizeof(unsigned char));
    memset(ctx->useq, 0, PDG_SEQ_CAP / 8 * sizeof(unsigned char));

    /* Allocate node and gene arrays */
    ctx->nodes = (struct _node *)pdg_alloc(ctx, STT_NOD * sizeof(struct _node));
    ctx->node_cap = STT_NOD;
    ctx->genes = (struct _gene *)pdg_alloc(ctx, MAX_GENES * sizeof(struct _gene));
    ctx->gene_cap = MAX_GENES;
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
    pdg_free(ctx, ctx->masks);

    free_meta(ctx);

    free(ctx);  /* Context itself always uses system malloc */
}

/*******************************************************************************
    Working arrays
*******************************************************************************/

/* A zeroed array of n elements of size bytes, holding the first keep of
   old's elements, in place of old (freed). NULL if out of memory, and old
   is kept. */
static void *grow(prodigal_ctx_t *ctx, void *old, size_t n, size_t size,
                  size_t keep) {
    void *p = n > SIZE_MAX / size ? NULL : pdg_calloc(ctx, n * size);
    if (p == NULL) return NULL;
    if (keep > 0) memcpy(p, old, keep * size);
    pdg_free(ctx, old);
    return p;
}

/*******************************************************************************
    Sequence encoding, shared with the CLI's readers (prodigal_internal.h)
*******************************************************************************/

/* Encode base ch at position pos (A 00, G 01, C 10, T 11; any other letter
   is ambiguous: C, flagged in useq). Returns 1 for G or C. */
static int encode_base(unsigned char *seq, unsigned char *useq, int pos,
                       char ch) {
    if (ch == 'g' || ch == 'G') { set(seq, 2 * pos); return 1; }
    if (ch == 't' || ch == 'T') { set(seq, 2 * pos); set(seq, 2 * pos + 1); }
    else if (ch == 'c' || ch == 'C') { set(seq, 2 * pos + 1); return 1; }
    else if (ch != 'a' && ch != 'A') { set(seq, 2 * pos + 1); set(useq, pos); }
    return 0;
}

/* -m: before position pos gets base ch, record the run of Ns that ch ends,
   if at least MASK_SIZE long. Returns 0, or -1 if out of memory. */
static int track_mask(prodigal_ctx_t *ctx, int *mask_beg, int pos, char ch) {
    int is_n = ch == 'N' || ch == 'n';
    mask *m;
    if (*mask_beg != -1 && !is_n) {
        if (pos - *mask_beg >= MASK_SIZE) {
            if (ctx->nmask == ctx->mask_cap) {
                m = (mask *)grow(ctx, ctx->masks,
                                 ctx->mask_cap > 0 ? 2 * (size_t)ctx->mask_cap : 64,
                                 sizeof(mask), (size_t)ctx->nmask);
                if (m == NULL) return -1;
                ctx->masks = m;
                ctx->mask_cap = ctx->mask_cap > 0 ? 2 * ctx->mask_cap : 64;
            }
            ctx->masks[ctx->nmask].begin = *mask_beg;
            ctx->masks[ctx->nmask].end = pos - 1;
            ctx->nmask++;
        }
        *mask_beg = -1;
    }
    if (*mask_beg == -1 && is_n) *mask_beg = pos;
    return 0;
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

void pdg_encode_begin(prodigal_ctx_t *ctx, pdg_encoder_t *e, int limit) {
    pdg_clear_sequence(ctx);
    e->len = 0;
    e->gc = 0;
    e->mask_beg = -1;
    e->limit = limit;
}

int pdg_encode_chars(prodigal_ctx_t *ctx, pdg_encoder_t *e, const char *s,
                     size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        if (s[i] < 'A' || s[i] > 'z') continue;
        if (e->len >= e->limit) return PRODIGAL_ERR_INVALID_INPUT;
        if (ctx->config.mask_regions &&
            track_mask(ctx, &e->mask_beg, e->len, s[i]) != 0)
            return PRODIGAL_ERR_NOMEM;
        e->gc += encode_base(ctx->seq, ctx->useq, e->len++, s[i]);
    }
    return PRODIGAL_OK;
}

int pdg_encode_spacer(prodigal_ctx_t *ctx, pdg_encoder_t *e) {
    int i;
    if (e->len + 12 > e->limit) return PRODIGAL_ERR_INVALID_INPUT;
    for (i = 0; i < 12; i++, e->len++)
        if (i % 4 < 2) encode_base(ctx->seq, ctx->useq, e->len, 'T');
    return PRODIGAL_OK;
}

void pdg_encode_end(prodigal_ctx_t *ctx, pdg_encoder_t *e) {
    ctx->slen = e->len;
    ctx->gc = e->len > 0 ? (double)e->gc / (double)e->len : 0.0;
    rcom_seq(ctx->seq, ctx->rseq, ctx->useq, e->len);
}

/*******************************************************************************
    Sequence input
*******************************************************************************/

/* Encode seqs[0, n) (NULL or empty ones add nothing), with a spacer after
   each if n > 1, as the sequence to work on; header names it. */
static int encode_seqs(prodigal_ctx_t *ctx, const char *const *seqs,
                       const int32_t *lens, int32_t n, const char *header) {
    pdg_encoder_t e;
    int32_t s;
    int rc = PRODIGAL_OK;

    pdg_encode_begin(ctx, &e, MAX_BASES);
    for (s = 0; s < n && rc == PRODIGAL_OK; s++) {
        if (seqs[s] != NULL && lens[s] > 0)
            rc = pdg_encode_chars(ctx, &e, seqs[s], (size_t)lens[s]);
        if (rc == PRODIGAL_OK && n > 1) rc = pdg_encode_spacer(ctx, &e);
    }
    if (rc != PRODIGAL_OK || e.len == 0) {
        ctx->slen = e.len;   /* clear what was encoded */
        pdg_clear_sequence(ctx);
        if (rc == PRODIGAL_ERR_NOMEM)
            PDG_SET_ERROR(ctx, rc, "Out of memory for masked regions");
        else if (rc != PRODIGAL_OK)
            PDG_SET_ERROR(ctx, rc, "Sequence too long (at most %d bases, not "
                          "counting characters outside A-z)", MAX_BASES);
        else
            PDG_SET_ERROR(ctx, rc = PRODIGAL_ERR_INVALID_INPUT,
                          "No valid bases found");
        return rc;
    }
    pdg_encode_end(ctx, &e);
    if (header != NULL) {
        strncpy(ctx->cur_header, header, MAX_LINE - 1);
        ctx->cur_header[MAX_LINE - 1] = '\0';
    }
    return PRODIGAL_OK;
}

int prodigal_set_sequence(prodigal_ctx_t *ctx, const char *seq, int32_t len,
                          const char *header) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (seq == NULL || len <= 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Sequence is NULL or length <= 0");
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    return encode_seqs(ctx, &seq, &len, 1, header);
}

int prodigal_set_training_sequences(prodigal_ctx_t *ctx, const char **seqs,
                                    const char **headers, const int32_t *lens,
                                    int32_t n_seqs) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (seqs == NULL || lens == NULL || n_seqs <= 0) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Invalid training sequence arguments");
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    return encode_seqs(ctx, seqs, lens, n_seqs,
                       headers != NULL ? headers[0] : NULL);
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

/* Room for at least cap nodes. A new array is zeroed and replaces the old
   one, whose contents are dropped. Returns 0, or -1 if out of memory (the
   caller sets the error). */
static int reserve_nodes(prodigal_ctx_t *ctx, int cap) {
    struct _node *nodes;
    if (cap <= ctx->node_cap) return 0;
    nodes = (struct _node *)grow(ctx, ctx->nodes, (size_t)cap,
                                 sizeof(struct _node), 0);
    if (nodes == NULL) return -1;
    ctx->nodes = nodes;
    ctx->nn = 0;
    ctx->node_cap = cap;
    return 0;
}

/* The sorted nodes of the current sequence under model tinf, in ctx->nodes.
   add_nodes and the steps after it expect zeroed nodes (add_nodes sets edge
   only when it is 1, and a stale elim drops a gene), so the previous list,
   ctx->nodes[0, ctx->nn), is cleared first; everything past it is zero.
   The array is first sized for slen/8 nodes, upstream's budget, which real
   genomes stay well under; low-complexity input can need up to 2*slen+2
   (add_nodes adds at most one node per pass), so when the nodes don't fit
   the array doubles and add_nodes runs again. */
static int build_nodes(prodigal_ctx_t *ctx, struct _training *tinf) {
    int nn;
    if (reserve_nodes(ctx, ctx->slen / 8) != 0) {
        /* Upstream's message, which the CLI prints */
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM, "Realloc failed on nodes");
        return PRODIGAL_ERR_NOMEM;
    }
    memset(ctx->nodes, 0, ctx->nn * sizeof(struct _node));
    while ((nn = add_nodes(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                           ctx->node_cap, ctx->config.closed_ends,
                           ctx->masks, ctx->nmask, tinf)) < 0) {
        if (reserve_nodes(ctx, 2 * ctx->node_cap) != 0) {
            /* add_nodes filled the array: clear it for the next call */
            memset(ctx->nodes, 0, (size_t)ctx->node_cap * sizeof(struct _node));
            ctx->nn = 0;
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                          "Out of memory for more than %d nodes in a %d bp "
                          "sequence (low-complexity input has up to 2 per bp)",
                          ctx->node_cap, ctx->slen);
            return PRODIGAL_ERR_NOMEM;
        }
    }
    ctx->nn = nn;
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
    gc_frame = calc_most_gc_frame(ctx, ctx->seq, ctx->slen);
    if (gc_frame == NULL) {
        /* Upstream's message, which the CLI prints */
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM, "Malloc failed on gc frame plot");
        return PRODIGAL_ERR_NOMEM;
    }
    record_gc_bias(gc_frame, ctx->nodes, ctx->nn, &ctx->tinf);
    pdg_free(ctx, gc_frame);
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
    if (ctx->tinf.uses_sd == 0) {
        struct _motif_counts *mc = (struct _motif_counts *)
            pdg_alloc(ctx, sizeof(struct _motif_counts));
        if (mc == NULL) {
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                          "Out of memory for the non-SD start motif counts");
            return PRODIGAL_ERR_NOMEM;
        }
        train_starts_nonsd(ctx->seq, ctx->rseq, ctx->slen, ctx->nodes,
                           ctx->nn, &ctx->tinf, mc);
        pdg_free(ctx, mc);
    }
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

/* Also prodigal_load_training's check. Only the fields training sets from
   the config and counts are checked: other doubles can be NaN
   (upstream-behavior.md, quirk 11). */
int prodigal_training_info(const void *data, size_t len,
                           prodigal_training_info_t *info, char *msg,
                           size_t msg_len) {
    char buf[256] = "";
    const char *why;
    int table = 0, uses_sd = 0;
    double gc = 0.0, st_wt = 0.0;
    if (info == NULL) {
        snprintf(buf, sizeof(buf), "info is NULL");
    } else if (data == NULL) {
        snprintf(buf, sizeof(buf), "Training data is NULL");
    } else if (len != sizeof(struct _training)) {
        snprintf(buf, sizeof(buf), "Invalid training data (expected %zu bytes, "
                 "got %zu)", sizeof(struct _training), len);
    } else {
#define FIELD(f, v) memcpy(&(v), (const char *)data + \
                           offsetof(struct _training, f), sizeof(v))
        FIELD(trans_table, table);
        FIELD(gc, gc);
        FIELD(st_wt, st_wt);
        FIELD(uses_sd, uses_sd);
#undef FIELD
        if ((why = pdg_trans_table_unsupported(table)) != NULL)
            snprintf(buf, sizeof(buf), "Training data has translation table "
                     "%d, which" PDG_TABLE_UNSUPPORTED, table, why);
        else if (!is_valid_gc(gc))
            snprintf(buf, sizeof(buf), "Training data has GC %g, outside [0, 1]",
                     gc);
        else if (!is_valid_start_weight(st_wt))
            snprintf(buf, sizeof(buf), "Training data has start weight %g, not "
                     "finite and > 0", st_wt);
        else if (uses_sd != 0 && uses_sd != 1)
            snprintf(buf, sizeof(buf), "Training data has SD flag %d, not 0 or 1",
                     uses_sd);
    }
    if (msg != NULL && msg_len > 0) snprintf(msg, msg_len, "%s", buf);
    if (buf[0] != '\0') return PRODIGAL_ERR_INVALID_INPUT;
    info->translation_table = table;
    info->uses_sd = uses_sd;
    info->gc_content = gc;
    info->start_weight = st_wt;
    return PRODIGAL_OK;
}

int prodigal_load_training(prodigal_ctx_t *ctx, const void *data, size_t len) {
    prodigal_training_info_t info;
    char why[256];
    int rc;
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    /* Checked before copying, so a rejected blob changes nothing */
    rc = prodigal_training_info(data, len, &info, why, sizeof(why));
    if (rc != PRODIGAL_OK) {
        PDG_SET_ERROR(ctx, rc, "%s", why);
        return rc;
    }
    memcpy(&ctx->tinf, data, sizeof(struct _training));
    ctx->trained = 1;
    return PRODIGAL_OK;
}

int prodigal_export_training(const prodigal_ctx_t *ctx, void **data_out,
                             size_t *len_out) {
    void *buf;
    if (ctx == NULL || data_out == NULL || len_out == NULL || !ctx->trained)
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
    const char *why;
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if ((why = pdg_trans_table_unsupported(table)) != NULL) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Translation table %d" PDG_TABLE_UNSUPPORTED, table, why);
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    ctx->config.trans_table = ctx->tinf.trans_table = table;
    return PRODIGAL_OK;
}

int prodigal_set_start_weight(prodigal_ctx_t *ctx, double weight) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (!is_valid_start_weight(weight)) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "Start weight must be finite and > 0 (got %g)", weight);
        return PRODIGAL_ERR_INVALID_INPUT;
    }
    ctx->config.start_weight = ctx->tinf.st_wt = weight;
    return PRODIGAL_OK;
}

int prodigal_set_gc(prodigal_ctx_t *ctx, double gc) {
    if (ctx == NULL) return PRODIGAL_ERR_INVALID_INPUT;
    if (!is_valid_gc(gc)) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INVALID_INPUT,
                      "GC must be in [0, 1] (got %g)", gc);
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
    struct _gene *genes;
    int ng;
    eliminate_bad_genes(ctx->nodes, ipath, tinf);
    /* When the genes don't fit (upstream truncated at MAX_GENES), the array
       doubles and add_genes runs again */
    while ((ng = add_genes(ctx->genes, ctx->gene_cap, ctx->nodes, ipath)) < 0) {
        genes = (struct _gene *)grow(ctx, ctx->genes, 2 * (size_t)ctx->gene_cap,
                                     sizeof(struct _gene), 0);
        if (genes == NULL) {
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                          "Out of memory for more than %d genes", ctx->gene_cap);
            return PRODIGAL_ERR_NOMEM;
        }
        ctx->genes = genes;
        ctx->gene_cap *= 2;
    }
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
            rc = keep_genes(ctx, ipath, tinf, seqnum);
            if (rc < 0) return rc;
            res->n_genes = rc;
        }
    }
    return PRODIGAL_OK;
}

/*******************************************************************************
    SOA allocation and extraction
*******************************************************************************/

#define ALIGN16(x) (((x) + 15) & ~(size_t)15)

/* The output for ctx->genes[0, ng), as record_gene_data captured them,
   with their proteins under model tinf if the config wants them. Sets
   *out, or the error. */
static int extract_soa(prodigal_ctx_t *ctx, int ng, struct _training *tinf,
                       prodigal_genes_soa_t **out) {
    const struct _gene *genes = ctx->genes;
    prodigal_genes_soa_t *soa;
    size_t n = (size_t)(ng > 0 ? ng : 1);  /* at least 1 to avoid zero alloc */
    size_t offset = 0, letters = 0;
    int proteins = ctx->config.want_proteins != 0;
    char *base;
    struct _rbs_text *text;
    int i, k;

    /* Unreachable (about MAX_BASES / 3 letters; genes overlap little): it
       guards the int32 offsets */
    for (i = 0; proteins && i < ng; i++)
        letters += (size_t)pdg_protein_length(&genes[i]);
    if (letters > INT32_MAX) {
        PDG_SET_ERROR(ctx, PRODIGAL_ERR_INTERNAL,
                      "%zu protein letters overflow 32-bit offsets", letters);
        return PRODIGAL_ERR_INTERNAL;
    }

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
    size_t sz_offsets  = proteins ? ALIGN16((n + 1) * sizeof(int32_t)) : 0;
    size_t sz_letters  = ALIGN16(letters);  /* aligned_alloc: whole blocks */

    size_t total = sz_begin + sz_end + sz_strand + sz_pleft + sz_pright +
                   sz_stype + sz_cscore + sz_sscore + sz_rscore + sz_uscore +
                   sz_tscore + sz_conf + sz_gc + sz_motif + sz_spacer +
                   sz_text + sz_offsets + sz_letters;

    soa = (prodigal_genes_soa_t *)malloc(sizeof(prodigal_genes_soa_t));
    if (soa == NULL) goto nomem;
    memset(soa, 0, sizeof(*soa));
    soa->n_genes = ng;

    if (ng == 0) {
        soa->_base = NULL;
        *out = soa;
        return PRODIGAL_OK;
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
        goto nomem;
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
    text              = (struct _rbs_text *)(base + offset); offset += sz_text;
    if (proteins) {
        soa->protein_offsets = (int32_t *)(base + offset); offset += sz_offsets;
        soa->protein_data    = base + offset;
        soa->protein_offsets[0] = 0;
    }

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
        if (proteins) {
            char *p = soa->protein_data + soa->protein_offsets[i];
            int len = pdg_protein_length(g);
            for (k = 0; k < len; k++)
                p[k] = pdg_protein_letter(g, k, ctx->seq, ctx->rseq, ctx->useq,
                                          ctx->slen, tinf);
            soa->protein_offsets[i + 1] = soa->protein_offsets[i] + len;
        }
    }

    *out = soa;
    return PRODIGAL_OK;

nomem:
    PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM, "Failed to allocate SOA output");
    return PRODIGAL_ERR_NOMEM;
}

/*******************************************************************************
    Gene finding — public API
*******************************************************************************/

/* tinf is the model that called the genes: bin's in metagenomic mode, the
   context's in single mode (bin -1). */
static void fill_stats(const prodigal_ctx_t *ctx, int ng, int nn, int bin,
                       const struct _training *tinf, prodigal_stats_t *stats) {
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
    struct _training *tinf;
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
        if (ng < 0) return ng;
        nn = ctx->nn;
    }

    tinf = (bin >= 0) ? ctx->meta[bin].tinf : &ctx->tinf;
    rc = extract_soa(ctx, ng, tinf, genes_out);
    if (rc != PRODIGAL_OK) return rc;
    fill_stats(ctx, ng, nn, bin, tinf, stats_out);
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
        char *prot;
        size_t n = (size_t)soa->n_genes, letters = 0;
        if (soa->protein_offsets != NULL)
            letters = (size_t)soa->protein_offsets[n] + n;   /* with NULs */
        /* The genes, their motif and spacer text, then their proteins */
        aos->_base = malloc(n * (sizeof(prodigal_gene_t) +
                                 sizeof(struct _rbs_text)) + letters);
        if (aos->_base == NULL) {
            free(aos);
            prodigal_genes_free(soa);
            PDG_SET_ERROR(ctx, PRODIGAL_ERR_NOMEM,
                          "Failed to allocate AOS output");
            return PRODIGAL_ERR_NOMEM;
        }
        aos->genes = (prodigal_gene_t *)aos->_base;
        text = (struct _rbs_text *)(aos->genes + n);
        prot = letters ? (char *)(text + n) : NULL;

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
            aos->genes[i].protein = prot;
            if (prot != NULL) {
                size_t len = (size_t)(soa->protein_offsets[i + 1] -
                                      soa->protein_offsets[i]);
                memcpy(prot, soa->protein_data + soa->protein_offsets[i], len);
                prot[len] = '\0';
                prot += len + 1;
            }
        }
    }

    prodigal_genes_free(soa);
    *genes_out = aos;
    return PRODIGAL_OK;
}

/*******************************************************************************
    ABI checks
*******************************************************************************/

int prodigal_abi_version(void) { return PRODIGAL_ABI_VERSION; }

/* The fields of each public struct, in declaration order: F(T, field) for
   a scalar or pointer, A(T, field) for an array. abi_layout is made from
   these lists, and so is a positional initializer of each struct
   (abi_complete_*), which gcc and clang reject when a list misses a field
   (-Wmissing-field-initializers, an error there); test_abi_layout checks
   the order and that every public struct is here. */
#define CONFIG_FIELDS(F, A, T) F(T, struct_size) F(T, trans_table) \
    F(T, closed_ends) F(T, mask_regions) F(T, force_nonsd) F(T, meta_mode) \
    F(T, start_weight) F(T, alloc_fn) F(T, free_fn) \
    F(T, allocator_user_data) F(T, log_callback) F(T, log_user_data) \
    F(T, progress_callback) F(T, progress_user_data) F(T, want_proteins)
#define SEQ_INFO_FIELDS(F, A, T) F(T, length) F(T, gc_content)
#define SOA_FIELDS(F, A, T) F(T, n_genes) F(T, begin) F(T, end) \
    F(T, strand) F(T, partial_left) F(T, partial_right) F(T, start_type) \
    F(T, cscore) F(T, sscore) F(T, rscore) F(T, uscore) F(T, tscore) \
    F(T, confidence) F(T, gc_cont) F(T, rbs_motif) F(T, rbs_spacer) \
    F(T, protein_offsets) F(T, protein_data) F(T, _base)
#define GENE_FIELDS(F, A, T) F(T, begin) F(T, end) F(T, strand) \
    F(T, partial_left) F(T, partial_right) F(T, start_type) F(T, cscore) \
    F(T, sscore) F(T, rscore) F(T, uscore) F(T, tscore) F(T, confidence) \
    F(T, gc_cont) F(T, rbs_motif) F(T, rbs_spacer) F(T, protein)
#define GENES_FIELDS(F, A, T) F(T, n_genes) F(T, genes) F(T, _base)
#define STATS_FIELDS(F, A, T) F(T, n_genes) F(T, n_nodes) F(T, gc_content) \
    F(T, translation_table) F(T, uses_sd) F(T, best_meta_bin) \
    A(T, best_meta_desc)
#define TRAINING_INFO_FIELDS(F, A, T) F(T, translation_table) \
    F(T, uses_sd) F(T, gc_content) F(T, start_weight)
#define ABI_FIELD_FIELDS(F, A, T) F(T, type) F(T, field) F(T, offset) \
    F(T, size)
#define PDG_ABI_STRUCTS(S) \
    S(prodigal_config_t, CONFIG_FIELDS) \
    S(prodigal_seq_info_t, SEQ_INFO_FIELDS) \
    S(prodigal_genes_soa_t, SOA_FIELDS) \
    S(prodigal_gene_t, GENE_FIELDS) \
    S(prodigal_genes_t, GENES_FIELDS) \
    S(prodigal_stats_t, STATS_FIELDS) \
    S(prodigal_training_info_t, TRAINING_INFO_FIELDS) \
    S(prodigal_abi_field_t, ABI_FIELD_FIELDS)

#define ROW(T, f) { #T, #f, offsetof(T, f), sizeof(((T *)0)->f) },
#define STRUCT_ROWS(T, FIELDS) { #T, NULL, 0, sizeof(T) }, FIELDS(ROW, ROW, T)
static const prodigal_abi_field_t abi_layout[] = {
    PDG_ABI_STRUCTS(STRUCT_ROWS)
    { NULL, NULL, 0, 0 }
};
#undef STRUCT_ROWS
#undef ROW

/* Each struct initialized from its list, in a type that is never used:
   the declarations only have to compile. (Like every universal zero
   initializer, a list of one field would not be checked; none has one.) */
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wmissing-field-initializers"
#endif
#define ZERO(T, f) 0,
#define ZERO_ARRAY(T, f) {0},
#define COMPLETE(T, FIELDS) \
    typedef char abi_complete_##T[sizeof((T){ FIELDS(ZERO, ZERO_ARRAY, T) })];
PDG_ABI_STRUCTS(COMPLETE)
#undef COMPLETE
#undef ZERO_ARRAY
#undef ZERO
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

const prodigal_abi_field_t *prodigal_abi_layout(size_t *n) {
    if (n != NULL) *n = sizeof(abi_layout) / sizeof(abi_layout[0]) - 1;
    return abi_layout;
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
