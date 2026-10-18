/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Internal header -- not part of the public API.
*******************************************************************************/

#ifndef PRODIGAL_INTERNAL_H
#define PRODIGAL_INTERNAL_H

#include <stdlib.h>
#include <string.h>
#include "prodigal.h"
#include "sequence.h"
#include "node.h"
#include "gene.h"
#include "training.h"
#include "metagenomic.h"
#include "dprog.h"
#include "bitmap.h"

/*******************************************************************************
    Internal context structure (opaque to public API consumers)
*******************************************************************************/

/* The most bases a sequence may have: the CLI's limit (it stops once the
   length plus MAX_LINE reaches MAX_SEQ). */
#define MAX_BASES (MAX_SEQ - MAX_LINE - 1)

/* The bases seq, rseq and useq hold: MAX_SEQ, and room for the spacer the
   CLI adds after training input it cut at the limit (upstream wrote it past
   its buffers). */
#define PDG_SEQ_CAP (MAX_SEQ + 16)

struct prodigal_ctx {
    prodigal_config_t config;         /* Snapshot of caller config */

    /* Training data */
    struct _training tinf;            /* ~558KB training model */
    int trained;                      /* Nonzero if training is loaded/complete */

    /* Sequence buffers */
    unsigned char *seq;               /* 2-bit encoded forward sequence */
    unsigned char *rseq;              /* 2-bit encoded reverse complement */
    unsigned char *useq;              /* Ambiguity bitmap (N bases) */
    int slen;                         /* Current sequence length in bp */
    int node_cap;                     /* Nodes the nodes array holds */
    double gc;                        /* GC content of current sequence */
    char cur_header[MAX_LINE];        /* Current sequence header */

    /* Masking (-m): the regions, in sequence order; the list grows as
       needed (upstream stopped at 5,000) */
    mask *masks;
    int nmask;
    int mask_cap;

    /* Node and gene working arrays */
    struct _node *nodes;
    int nn;                           /* Current node count */
    struct _gene *genes;
    int gene_cap;                     /* Genes the genes array holds */

    /* Metagenomic bins (lazily initialized) */
    struct _metagenomic_bin *meta;

    /* Error state */
    char error_msg[1024];
    int error_code;
};

/*******************************************************************************
    Internal helpers
*******************************************************************************/

static inline void *pdg_alloc(prodigal_ctx_t *ctx, size_t size) {
    if (ctx->config.alloc_fn)
        return ctx->config.alloc_fn(size, ctx->config.allocator_user_data);
    return malloc(size);
}

/* Zeroed memory; from calloc without a custom allocator, so large blocks
   are committed only as they are touched. */
static inline void *pdg_calloc(prodigal_ctx_t *ctx, size_t size) {
    void *p;
    if (ctx->config.alloc_fn == NULL) return calloc(1, size);
    p = ctx->config.alloc_fn(size, ctx->config.allocator_user_data);
    if (p != NULL) memset(p, 0, size);
    return p;
}

static inline void pdg_free(prodigal_ctx_t *ctx, void *ptr) {
    if (ptr == NULL) return;
    if (ctx->config.free_fn)
        ctx->config.free_fn(ptr, ctx->config.allocator_user_data);
    else
        free(ptr);
}

/* The format is the first of the variadic arguments, so a message without
   arguments is standard C99 (no GNU ##__VA_ARGS__). */
#define PDG_LOG(ctx, ...) do { \
    if ((ctx)->config.log_callback) { \
        char _pdg_buf[512]; \
        snprintf(_pdg_buf, sizeof(_pdg_buf), __VA_ARGS__); \
        (ctx)->config.log_callback(_pdg_buf, (ctx)->config.log_user_data); \
    } \
} while(0)

#define PDG_SET_ERROR(ctx, code, ...) do { \
    (ctx)->error_code = (code); \
    snprintf((ctx)->error_msg, sizeof((ctx)->error_msg), __VA_ARGS__); \
} while(0)

/*******************************************************************************
    The pipeline, shared by the CLI (main.c) and the API (prodigal_api.c), so
    the two compute the same genes (docs/internals.md, "One pipeline"). Each
    function returns PRODIGAL_OK or an error code with the context's error
    message set, unless noted.
*******************************************************************************/

/* Clear the sequence state (buffers, length, masks, GC, header and error)
   for the next sequence. The nodes are build_nodes's (prodigal_api.c). */
void pdg_clear_sequence(prodigal_ctx_t *ctx);

/* Encoding a sequence into ctx (seq, useq, and with -m the masks), shared
   by the library's prodigal_set_sequence / prodigal_set_training_sequences
   and the CLI's readers (sequence.c), so both keep, drop and mask the same
   characters. pdg_encode_begin clears the sequence; pdg_encode_chars
   encodes the characters in 'A'..'z' of s[0, n) and drops the rest;
   pdg_encode_spacer adds TTAATTAATTAA (a stop in every frame) between
   training sequences; pdg_encode_end sets slen and gc and computes the
   reverse strand. The two can fail with PRODIGAL_ERR_INVALID_INPUT past
   limit bases, or PRODIGAL_ERR_NOMEM if the mask list can't grow; they
   set no message. */
typedef struct {
    int len;                          /* bases encoded */
    int gc;                           /* of them G or C */
    int mask_beg;                     /* -m: start of the current N run, or -1 */
    int limit;                        /* the most bases allowed */
} pdg_encoder_t;
void pdg_encode_begin(prodigal_ctx_t *ctx, pdg_encoder_t *e, int limit);
int pdg_encode_chars(prodigal_ctx_t *ctx, pdg_encoder_t *e, const char *s,
                     size_t n);
int pdg_encode_spacer(prodigal_ctx_t *ctx, pdg_encoder_t *e);
void pdg_encode_end(prodigal_ctx_t *ctx, pdg_encoder_t *e);

/* Training on the sequence in ctx, in stages. report, if not NULL, is
   called before (done 0) and after (done 1) each stage: prodigal_train logs
   from it, and the CLI prints upstream's progress text. */
enum { PDG_TRAIN_NODES, PDG_TRAIN_GC_BIAS, PDG_TRAIN_INITIAL_GENES,
       PDG_TRAIN_CODING, PDG_TRAIN_STARTS };
typedef void (*pdg_train_report_fn)(const prodigal_ctx_t *ctx, int stage,
                                    int done);
int pdg_train(prodigal_ctx_t *ctx, pdg_train_report_fn report);

/* The nodes of the sequence under model tinf, built and scored. In single
   genome mode this is the first step (the CLI writes its start file, -s,
   after it), with ctx->tinf. In metagenomic mode the CLI calls it after
   pdg_meta with the best bin's model: GFF column 6 and the start file come
   from those nodes (docs/upstream-behavior.md, quirk 2). */
int pdg_score_nodes(prodigal_ctx_t *ctx, struct _training *tinf);

/* Single genome mode, after pdg_score_nodes: the genes, in ctx->genes, and
   their number, or PRODIGAL_ERR_NOMEM if the gene array can't grow to hold
   them. seqnum is the 1-based sequence number in the genes' IDs. */
int pdg_single_genes(prodigal_ctx_t *ctx, int seqnum);

/* Metagenomic mode: the bins, allocated on first use. */
int pdg_meta_init(prodigal_ctx_t *ctx);

/* Metagenomic mode: run every bin whose GC window holds the sequence (GC
   ctx->gc) and keep the genes of the best one in ctx->genes. Afterwards
   ctx->nodes holds the last bin's nodes, not the best bin's. */
typedef struct {
    int n_genes;   /* genes in ctx->genes */
    int bin;       /* index of the bin that called them */
    int n_nodes;   /* nodes under that bin's model */
} pdg_meta_result_t;
int pdg_meta(prodigal_ctx_t *ctx, int seqnum, pdg_meta_result_t *res);

#endif /* PRODIGAL_INTERNAL_H */
