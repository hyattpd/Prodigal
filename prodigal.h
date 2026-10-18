/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

/**
 * @file prodigal.h
 * @brief Reentrant C library API for Prodigal gene prediction.
 *
 * Overview
 * --------
 * Prodigal (PROkaryotic DynamIc Programming Genefinding ALgorithm) predicts
 * protein-coding genes in prokaryotic genomes. This header defines a library
 * API that allows embedding Prodigal in other applications without shelling
 * out to the command-line tool.
 *
 * Two usage patterns are supported:
 *
 *   High-level (single-call):
 *     config_init -> create -> set_sequence -> find_genes -> free -> destroy
 *
 *   Low-level (multi-step):
 *     config_init -> create -> set_training_sequences -> train ->
 *     set_sequence -> find_genes -> free -> (repeat per sequence) -> destroy
 *
 * Lifecycle
 * ---------
 *   1. Call prodigal_config_init() to populate a config with safe defaults.
 *   2. Modify config fields as needed (trans_table, meta_mode, callbacks...).
 *   3. Call prodigal_create() to allocate an opaque context.
 *   4. Load input:
 *      - Metagenomic mode: call prodigal_set_sequence() per contig, then
 *        prodigal_find_genes().  Training data is built-in (50 pre-trained
 *        models, lazily initialized on first use).
 *      - Single-genome mode: call prodigal_set_training_sequences() with all
 *        contigs, then prodigal_train().  Alternatively, load a previously
 *        exported training file with prodigal_load_training().  Then call
 *        prodigal_set_sequence() + prodigal_find_genes() per contig.
 *   5. Free each output with prodigal_genes_free() / prodigal_genes_aos_free().
 *   6. Call prodigal_destroy() when done.
 *
 * Memory model
 * ------------
 *   - The context struct itself is allocated with system malloc/free.
 *   - Internal working buffers (sequence, nodes, genes) use the custom
 *     allocator if provided, otherwise system malloc/free.
 *   - Output structs (prodigal_genes_soa_t, prodigal_genes_t) are always
 *     allocated with system malloc so they can be freed without a context
 *     reference.  This matches the convention used by GPL-boundary's
 *     FastTree integration.
 *   - prodigal_destroy(NULL) and prodigal_genes_free(NULL) are safe no-ops.
 *
 * Output formats
 * --------------
 *   - SOA (prodigal_genes_soa_t): Structure-of-Arrays layout ideal for
 *     columnar / Apache Arrow consumers.  All numeric arrays are carved from
 *     a single 16-byte-aligned backing allocation (_base).  The strings
 *     (rbs_motif, rbs_spacer) live in it too: they are valid until
 *     prodigal_genes_free() and are not freed on their own.
 *   - AOS (prodigal_genes_t): Array-of-Structures layout for convenient
 *     per-gene iteration.  Also backed by a single allocation.
 *
 * Error handling
 * --------------
 *   - All fallible functions return int: PRODIGAL_OK (0) on success,
 *     negative error codes on failure.
 *   - prodigal_strerror() maps an error code to a static category string.
 *   - prodigal_last_error() returns a detailed message from the context,
 *     valid until the next API call on that context or prodigal_destroy().
 *   - The context is reusable after an error -- call prodigal_set_sequence()
 *     again and retry.
 *   - The library never calls exit(), abort(), or assert().
 *   - The library never writes to stdout or stderr.  Use log_callback for
 *     diagnostic output.
 *
 * Thread safety
 * -------------
 *   - Each context is independent.  Multiple contexts can be used
 *     concurrently from different threads without synchronization.
 *   - A single context must not be used from multiple threads simultaneously.
 *   - No global mutable state exists in the library.
 *
 * GPL-boundary integration
 * ------------------------
 *   - prodigal_config_t uses struct_size as its first field for ABI
 *     versioning.  The Rust adapter checks this at test time.
 *   - Build with -DPRODIGAL_NO_MAIN to exclude the CLI entry point when
 *     compiling as a static library for embedding.
 *   - The Makefile produces libprodigal.a (static) and libprodigal.so
 *     (shared).  The shared library has no zlib dependency; gzip support
 *     is only used by the CLI's file I/O path.
 *
 * Quick example (metagenomic mode)
 * --------------------------------
 * @code
 *   prodigal_config_t config;
 *   prodigal_config_init(&config);
 *   config.meta_mode = 1;
 *
 *   prodigal_ctx_t *ctx = prodigal_create(&config);
 *   if (ctx == NULL) { handle error }
 *
 *   prodigal_set_sequence(ctx, seq_chars, seq_len, "contig_1");
 *
 *   prodigal_genes_soa_t *genes = NULL;
 *   prodigal_stats_t stats;
 *   int rc = prodigal_find_genes(ctx, &genes, &stats);
 *   if (rc != PRODIGAL_OK) { handle error }
 *
 *   for (int i = 0; i < genes->n_genes; i++) {
 *       printf("gene %d: %d..%d strand=%d score=%.1f\n", i,
 *              genes->begin[i], genes->end[i], genes->strand[i],
 *              genes->cscore[i] + genes->sscore[i]);
 *   }
 *
 *   prodigal_genes_free(genes);
 *   prodigal_destroy(ctx);
 * @endcode
 *
 * Quick example (single-genome with training)
 * --------------------------------------------
 * @code
 *   prodigal_config_t config;
 *   prodigal_config_init(&config);
 *
 *   prodigal_ctx_t *ctx = prodigal_create(&config);
 *
 *   // Load all contigs for training
 *   prodigal_set_training_sequences(ctx, seqs, headers, lens, n_seqs);
 *   int rc = prodigal_train(ctx);
 *   if (rc != PRODIGAL_OK) { handle error }
 *
 *   // Find genes per contig
 *   for (int s = 0; s < n_seqs; s++) {
 *       prodigal_set_sequence(ctx, seqs[s], lens[s], headers[s]);
 *       prodigal_genes_soa_t *genes = NULL;
 *       prodigal_find_genes(ctx, &genes, NULL);
 *       // ... process genes ...
 *       prodigal_genes_free(genes);
 *   }
 *
 *   prodigal_destroy(ctx);
 * @endcode
 *
 * Quick example (custom allocator)
 * --------------------------------
 * @code
 *   static void *my_alloc(size_t size, void *ud) { return my_pool_alloc(ud, size); }
 *   static void  my_free(void *ptr, void *ud)    { my_pool_free(ud, ptr); }
 *
 *   prodigal_config_t config;
 *   prodigal_config_init(&config);
 *   config.alloc_fn = my_alloc;
 *   config.free_fn  = my_free;
 *   config.allocator_user_data = my_pool;
 *   // ... create, use, destroy as usual ...
 * @endcode
 */

#ifndef PRODIGAL_API_H
#define PRODIGAL_API_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Symbol visibility.
   For static linking (the default), PRODIGAL_API is empty.
   For building a shared library (DLL), define PRODIGAL_BUILDING_DLL.
   For consuming a shared library on Windows, define PRODIGAL_DLL. */
#if defined(_WIN32) || defined(__CYGWIN__)
  #if defined(PRODIGAL_BUILDING_DLL)
    #define PRODIGAL_API __declspec(dllexport)
  #elif defined(PRODIGAL_DLL)
    #define PRODIGAL_API __declspec(dllimport)
  #else
    #define PRODIGAL_API
  #endif
#elif defined(__GNUC__) && __GNUC__ >= 4
  #define PRODIGAL_API __attribute__((visibility("default")))
#else
  #define PRODIGAL_API
#endif

/** @name Compile-time version constants */
/**@{*/
#define PRODIGAL_VERSION_MAJOR 2
#define PRODIGAL_VERSION_MINOR 6
#define PRODIGAL_VERSION_PATCH 4
#define PRODIGAL_VERSION_STRING "2.6.4"
/**@}*/

/**
 * @name Error codes
 * Zero is success; all errors are negative.  Use prodigal_strerror() for a
 * human-readable category and prodigal_last_error() for a detailed message.
 */
/**@{*/
#define PRODIGAL_OK                   0
#define PRODIGAL_ERR_NOMEM           -1  /**< Memory allocation failed */
#define PRODIGAL_ERR_INVALID_CONFIG  -2  /**< Bad config (NULL, wrong struct_size, invalid field) */
#define PRODIGAL_ERR_INVALID_INPUT   -3  /**< Bad input (NULL sequence, no sequence loaded) */
#define PRODIGAL_ERR_INTERNAL        -4  /**< Internal error (should not happen) */
#define PRODIGAL_ERR_SEQ_TOO_SHORT   -5  /**< Sequence < 20000 bp for single-genome training */
#define PRODIGAL_ERR_CANCELLED       -6  /**< Cancelled by progress_callback returning nonzero */
/**@}*/

/**
 * Opaque context.  All computation state lives here.
 * Created by prodigal_create(), destroyed by prodigal_destroy().
 */
typedef struct prodigal_ctx prodigal_ctx_t;

/*******************************************************************************
    Configuration
*******************************************************************************/

/**
 * Configuration struct passed to prodigal_create().
 *
 * Initialize with prodigal_config_init() which sets struct_size and all
 * defaults.  Modify fields as needed, then pass to prodigal_create().
 * The context takes a snapshot of the config; further changes to the config
 * struct after create() have no effect.
 *
 * New fields are always appended at the end.  Never reorder or remove fields.
 * ABI versioning is via struct_size (must be the first field).
 */
typedef struct {
    /** Must be first field.  Set by prodigal_config_init().  Used for ABI
     *  version detection: prodigal_create() rejects configs where
     *  struct_size != sizeof(prodigal_config_t). */
    size_t struct_size;

    /** NCBI translation table number.  Default: 11 (Standard Microbial).
     *  Valid: 1-6, 9-16, 21-25.  Invalid values cause prodigal_create()
     *  to return NULL. */
    int    trans_table;

    /** Nonzero: closed ends -- do not allow genes to run off sequence edges.
     *  Default: 0 (allow edge genes). */
    int    closed_ends;

    /** Nonzero: treat runs of >= 50 N's as masked regions; don't build genes
     *  across them.  Default: 0. */
    int    mask_regions;

    /** Nonzero: bypass Shine-Dalgarno trainer and force a full upstream
     *  motif scan.  Default: 0 (auto-detect SD usage). */
    int    force_nonsd;

    /** Nonzero: metagenomic mode (use 50 pre-trained models).
     *  Zero: single-genome mode (requires training).  Default: 0. */
    int    meta_mode;

    /** Start score weight.  Affects the balance between coding potential and
     *  start signal strength.  Default: 4.35.  Rarely needs changing. */
    double start_weight;

    /** Custom allocator for internal working buffers.  If NULL, system
     *  malloc is used.  Must return memory suitable for any alignment
     *  (16-byte aligned recommended).  The context struct itself and all
     *  output structs always use system malloc regardless of this setting. */
    void  *(*alloc_fn)(size_t size, void *user_data);

    /** Custom deallocator paired with alloc_fn.  If NULL, system free is
     *  used.  Called with pointers previously returned by alloc_fn. */
    void   (*free_fn)(void *ptr, void *user_data);

    /** User data pointer passed as the second argument to both alloc_fn
     *  and free_fn.  Typically a pool or arena handle. */
    void   *allocator_user_data;

    /** Logging callback.  Receives diagnostic messages (progress, warnings).
     *  If NULL, messages are silently discarded.  The msg string is valid
     *  only for the duration of the callback. */
    void   (*log_callback)(const char *msg, void *user_data);

    /** User data pointer passed to log_callback. */
    void   *log_user_data;

    /** Progress callback for long-running operations (metagenomic scoring).
     *  Called with a stage name and fraction done [0.0, 1.0].
     *  Return 0 to continue, nonzero to cancel (prodigal_find_genes will
     *  return PRODIGAL_ERR_CANCELLED).  If NULL, no progress reporting. */
    int    (*progress_callback)(const char *stage, double frac_done,
                                void *user_data);

    /** User data pointer passed to progress_callback. */
    void   *progress_user_data;
} prodigal_config_t;

/*******************************************************************************
    Sequence info
*******************************************************************************/

/** Basic information about the currently loaded sequence. */
typedef struct {
    /** Encoded sequence length in bp: the input's characters in 'A'..'z'
     *  (the rest are dropped), plus, after prodigal_set_training_sequences
     *  with several sequences, 12 for each spacer (one after each). */
    int32_t length;
    double  gc_content;         /**< GC fraction [0, 1].  0.0 if no sequence loaded. */
} prodigal_seq_info_t;

/*******************************************************************************
    Output: Structure of Arrays (SOA)
*******************************************************************************/

/**
 * Gene prediction results in Structure-of-Arrays layout.
 *
 * Preferred output format for columnar / Apache Arrow consumers.  All numeric
 * arrays are carved from a single 16-byte-aligned backing allocation pointed
 * to by _base.  Call prodigal_genes_free() to release.
 *
 * Every value is the one the CLI prints in its GFF attributes (partial=,
 * start_type=, rbs_motif=, ..., tscore=): the values that chose the gene.
 * In metagenomic mode the CLI's GFF column 6 differs from cscore + sscore
 * (an upstream quirk, docs/upstream-behavior.md).
 *
 * The strings (rbs_motif, rbs_spacer) are stored in the same allocation:
 * they are valid until prodigal_genes_free() and are not freed on their own.
 *
 * If n_genes == 0, all array pointers may be NULL and _base is NULL.
 */
typedef struct {
    int32_t  n_genes;           /**< Number of predicted genes */

    int32_t *begin;             /**< 1-based left coordinate, length n_genes */
    int32_t *end;               /**< 1-based right coordinate, length n_genes */
    int32_t *strand;            /**< +1 forward, -1 reverse, length n_genes */

    int32_t *partial_left;      /**< 1 if gene runs off left edge, length n_genes */
    int32_t *partial_right;     /**< 1 if gene runs off right edge, length n_genes */
    int32_t *start_type;        /**< 0=ATG, 1=GTG, 2=TTG, 3=Edge, length n_genes */

    double  *cscore;            /**< Coding score (6-mer log-odds), length n_genes */
    double  *sscore;            /**< Start score (tscore+rscore+uscore), length n_genes */
    double  *rscore;            /**< RBS motif score, length n_genes */
    double  *uscore;            /**< Upstream composition score, length n_genes */
    double  *tscore;            /**< Start type score, length n_genes */
    double  *confidence;        /**< Confidence in [50, 100], length n_genes */
    double  *gc_cont;           /**< Per-gene GC content, length n_genes */

    const char **rbs_motif;     /**< RBS motif, as the CLI prints it (e.g. "AGGAG", "TAA", "None"), length n_genes */
    const char **rbs_spacer;    /**< RBS spacer, as the CLI prints it (e.g. "5-10bp", "8bp", "None"), length n_genes */

    void    *_base;             /**< Single backing allocation (16-byte aligned).  Internal. */
} prodigal_genes_soa_t;

/*******************************************************************************
    Output: Array of Structures (AOS)
*******************************************************************************/

/** Per-gene data in struct form.  See prodigal_genes_soa_t for field docs. */
typedef struct {
    int32_t begin;
    int32_t end;
    int32_t strand;
    int32_t partial_left;
    int32_t partial_right;
    int32_t start_type;
    double  cscore;
    double  sscore;
    double  rscore;
    double  uscore;
    double  tscore;
    double  confidence;
    double  gc_cont;
    const char *rbs_motif;      /**< Valid until prodigal_genes_aos_free() */
    const char *rbs_spacer;     /**< Valid until prodigal_genes_aos_free() */
} prodigal_gene_t;

/**
 * Gene prediction results in Array-of-Structures layout.
 * Call prodigal_genes_aos_free() to release.
 */
typedef struct {
    int32_t         n_genes;    /**< Number of predicted genes */
    prodigal_gene_t *genes;     /**< Array of n_genes entries */
    void            *_base;     /**< Single backing allocation.  Internal. */
} prodigal_genes_t;

/*******************************************************************************
    Statistics
*******************************************************************************/

/**
 * Computation statistics.  Pointer-free, safe to memcpy.
 * Passed by pointer to prodigal_find_genes(); may be NULL if not needed.
 */
typedef struct {
    int32_t n_genes;            /**< Number of genes found */
    int32_t n_nodes;            /**< Number of start/stop nodes evaluated */
    double  gc_content;         /**< Sequence GC content */
    int32_t translation_table;  /**< Translation table used */
    int32_t uses_sd;            /**< Nonzero if Shine-Dalgarno motifs used */
    int32_t best_meta_bin;      /**< Best metagenomic bin index (-1 if not meta) */
    char    best_meta_desc[512];/**< Description of best metagenomic bin */
} prodigal_stats_t;

/*******************************************************************************
    API Functions
*******************************************************************************/

/**
 * Initialize a config struct with safe defaults.
 *
 * Sets struct_size, trans_table=11, start_weight=4.35, all other fields to
 * zero/NULL.  A config initialized this way is ready for prodigal_create()
 * without further modification (defaults to single-genome mode, table 11).
 *
 * @param config  Pointer to config struct to initialize.  Must not be NULL.
 */
PRODIGAL_API void prodigal_config_init(prodigal_config_t *config);

/**
 * Create a new Prodigal context.
 *
 * Allocates all internal buffers.  The config is snapshotted; the caller may
 * modify or discard it after this call.
 *
 * Returns NULL on failure (bad config, OOM, invalid translation table).
 * When NULL is returned, the caller cannot distinguish the cause via
 * last_error (no context exists).  Recheck inputs manually if needed.
 *
 * @param config  Pointer to initialized config.  Must not be NULL.
 * @return New context, or NULL on failure.
 */
PRODIGAL_API prodigal_ctx_t *prodigal_create(const prodigal_config_t *config);

/**
 * Destroy a context and free all associated memory.
 *
 * Safe to call with NULL (no-op).  After this call, the context pointer is
 * invalid.  Any output structs (prodigal_genes_soa_t, prodigal_genes_t)
 * previously returned remain valid until explicitly freed.
 *
 * @param ctx  Context to destroy, or NULL.
 */
PRODIGAL_API void prodigal_destroy(prodigal_ctx_t *ctx);

/**
 * Load a single sequence for gene finding.
 *
 * Encodes the raw ASCII nucleotide string into Prodigal's internal 2-bit
 * representation, computes the reverse complement, and calculates GC content.
 * Accepts A/C/G/T/N in upper or lower case.  Characters outside 'A'..'z'
 * (digits, '-', '.', '*', whitespace, ...) are dropped, as the CLI drops
 * them, which shifts the coordinates of everything after them:
 * prodigal_get_seq_info() reports how many characters were kept, so a
 * caller can compare that with len.  Other characters in 'A'..'z' (IUPAC
 * codes, '[', '_', ...) are ambiguous bases, encoded as C with an
 * ambiguity flag.
 *
 * At most 31,989,999 characters may be kept (MAX_SEQ - MAX_LINE - 1), the
 * CLI's limit; a longer sequence is an error (PRODIGAL_ERR_INVALID_INPUT).
 *
 * Resets internal state from any previous sequence.  The context can be
 * reused across multiple set_sequence + find_genes cycles.
 *
 * @param ctx     Context.
 * @param seq     Raw nucleotide string (not null-terminated required; len used).
 * @param len     Number of characters in seq.  Must be > 0.
 * @param header  Sequence name/header.  May be NULL.
 * @return PRODIGAL_OK or error code.
 */
PRODIGAL_API int prodigal_set_sequence(prodigal_ctx_t *ctx,
                                       const char *seq, int32_t len,
                                       const char *header);

/**
 * Load multiple sequences concatenated for single-genome training.
 *
 * Concatenates all sequences with TTAATTAATTAA stop-codon spacers (forcing
 * stops in all 6 reading frames), matching Prodigal's training behavior.
 * After this call, use prodigal_train() to build the model.  Characters are
 * kept or dropped as in prodigal_set_sequence(), and the same limit applies
 * to the whole, spacers included: a longer input is an error, where the
 * CLI warns and trains on the start of it.
 *
 * @param ctx      Context.
 * @param seqs     Array of n_seqs raw nucleotide strings.
 * @param headers  Array of n_seqs header strings (may be NULL).
 * @param lens     Array of n_seqs sequence lengths.
 * @param n_seqs   Number of sequences.  Must be > 0.
 * @return PRODIGAL_OK or error code.
 */
PRODIGAL_API int prodigal_set_training_sequences(prodigal_ctx_t *ctx,
                                                  const char **seqs,
                                                  const char **headers,
                                                  const int32_t *lens,
                                                  int32_t n_seqs);

/**
 * Get information about the currently loaded sequence.
 *
 * Returns length=0 and gc_content=0.0 if no sequence is loaded.
 *
 * @param ctx   Context.
 * @param info  Output struct to populate.
 * @return PRODIGAL_OK or error code.
 */
PRODIGAL_API int prodigal_get_seq_info(const prodigal_ctx_t *ctx,
                                       prodigal_seq_info_t *info);

/**
 * Train a gene-finding model from the loaded training sequences.
 *
 * Requires prior call to prodigal_set_training_sequences() with >= 20000 bp
 * of concatenated sequence.  Runs the full Prodigal training pipeline:
 * node creation, GC bias analysis, initial dynamic programming, dicodon
 * statistics, RBS scoring, and start codon weight training.
 *
 * Training starts from a fresh model: the translation table and start
 * weight are the context's (the config's, or the last value given to
 * prodigal_set_translation_table() or prodigal_set_start_weight()), the GC
 * is the training sequence's, and nothing else carries over from a model
 * trained, loaded or adjusted before.  If training fails after checking
 * the sequence, the context has no model until it is trained again or a
 * model is loaded.
 *
 * For the same input and options, the model equals the training file the
 * CLI writes with -t, provided the input is within the limits listed under
 * "Input limits" in docs/internals.md.
 *
 * After training, call prodigal_set_sequence() + prodigal_find_genes()
 * for each contig.  Alternatively, export the model with
 * prodigal_export_training() for later reuse.
 *
 * Not needed in metagenomic mode (training is built-in).
 *
 * @param ctx  Context with training sequences loaded.
 * @return PRODIGAL_OK, PRODIGAL_ERR_SEQ_TOO_SHORT, or other error code.
 */
PRODIGAL_API int prodigal_train(prodigal_ctx_t *ctx);

/**
 * Load a pre-trained model from a binary blob.
 *
 * The blob must be exactly sizeof(struct _training) bytes, as produced by
 * prodigal_export_training() or the CLI's -t flag.  Binary format is
 * architecture-specific (not portable across endianness or struct padding).
 * A rejected blob leaves the context unchanged.
 *
 * @param ctx   Context.
 * @param data  Pointer to training data blob.
 * @param len   Size of blob in bytes.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_load_training(prodigal_ctx_t *ctx,
                                        const void *data, size_t len);

/**
 * Export the current training model as a binary blob.
 *
 * The caller receives a malloc'd buffer that must be freed with free().
 * The blob is binary-compatible with Prodigal's -t training file format.
 *
 * @param ctx       Context with training data (from train() or load()).
 * @param data_out  Receives pointer to malloc'd blob.
 * @param len_out   Receives size of blob in bytes.
 * @return PRODIGAL_OK or error code.
 */
PRODIGAL_API int prodigal_export_training(const prodigal_ctx_t *ctx,
                                          void **data_out, size_t *len_out);

/** @name Training parameter setters
 *  Fine-grained control over individual training parameters.  These modify
 *  the current model directly.  Useful for experimentation or for tweaking
 *  a loaded model.  The translation table and start weight also apply to
 *  later training; prodigal_train() replaces the GC and SD usage, so set
 *  those after training or loading.
 */
/**@{*/

/**
 * Set the NCBI translation table.
 * @param table  Valid: 1-6, 9-16, 21-25.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_set_translation_table(prodigal_ctx_t *ctx, int table);

/**
 * Set the start score weight.
 * @param weight  Must be > 0.  Default: 4.35.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_set_start_weight(prodigal_ctx_t *ctx, double weight);

/**
 * Set the GC content in the training model.
 * @param gc  Must be in [0.0, 1.0].
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_set_gc(prodigal_ctx_t *ctx, double gc);

/**
 * Set whether the model uses Shine-Dalgarno motifs.
 * @param uses_sd  Nonzero: use SD.  Zero: use upstream motif scan.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_set_uses_sd(prodigal_ctx_t *ctx, int uses_sd);
/**@}*/

/**
 * Find genes in the currently loaded sequence (SOA output).
 *
 * In single-genome mode, requires prior training (prodigal_train() or
 * prodigal_load_training()).  In metagenomic mode, training is built-in
 * and lazily initialized on first call (~27 MB for 50 models).
 *
 * On success, *genes_out is set to a newly allocated SOA struct.  The caller
 * must free it with prodigal_genes_free().  If stats_out is non-NULL, it is
 * populated with computation statistics.
 *
 * On error, *genes_out is set to NULL.
 *
 * @param ctx        Context with a sequence loaded.
 * @param genes_out  Receives pointer to SOA output (caller frees).
 * @param stats_out  Receives computation stats, or NULL to skip.
 * @return PRODIGAL_OK or error code.
 */
PRODIGAL_API int prodigal_find_genes(prodigal_ctx_t *ctx,
                                     prodigal_genes_soa_t **genes_out,
                                     prodigal_stats_t *stats_out);

/**
 * Find genes in the currently loaded sequence (AOS output).
 *
 * Same as prodigal_find_genes() but returns an Array-of-Structures layout.
 * Free with prodigal_genes_aos_free().
 *
 * @param ctx        Context with a sequence loaded.
 * @param genes_out  Receives pointer to AOS output (caller frees).
 * @param stats_out  Receives computation stats, or NULL to skip.
 * @return PRODIGAL_OK or error code.
 */
PRODIGAL_API int prodigal_find_genes_aos(prodigal_ctx_t *ctx,
                                         prodigal_genes_t **genes_out,
                                         prodigal_stats_t *stats_out);

/**
 * Free a SOA gene output struct.
 * Safe to call with NULL (no-op).  Always uses system free().
 * @param genes  SOA struct to free, or NULL.
 */
PRODIGAL_API void prodigal_genes_free(prodigal_genes_soa_t *genes);

/**
 * Free an AOS gene output struct.
 * Safe to call with NULL (no-op).  Always uses system free().
 * @param genes  AOS struct to free, or NULL.
 */
PRODIGAL_API void prodigal_genes_aos_free(prodigal_genes_t *genes);

/**
 * Get a static human-readable string for an error code.
 *
 * Thread-safe; returns a pointer to a string literal.
 *
 * @param error_code  An error code (PRODIGAL_OK, PRODIGAL_ERR_*, etc.).
 * @return Static string like "Success", "Out of memory", "Unknown error".
 */
PRODIGAL_API const char *prodigal_strerror(int error_code);

/**
 * Get a detailed error message from the last failed operation.
 *
 * Returns a pointer to an internal buffer in the context.  Valid until the
 * next API call on the same context, or until prodigal_destroy().
 * Returns "" on a fresh context with no errors.
 * Returns "NULL context" if ctx is NULL.
 *
 * @param ctx  Context, or NULL.
 * @return Error message string (never NULL).
 */
PRODIGAL_API const char *prodigal_last_error(const prodigal_ctx_t *ctx);

/**
 * Get the library version string at runtime.
 *
 * Useful for FFI consumers that cannot use compile-time macros.
 * Returns a pointer to a string literal (e.g., "2.6.3").
 *
 * @return Static version string.
 */
PRODIGAL_API const char *prodigal_version_string(void);

#ifdef __cplusplus
}
#endif

#endif /* PRODIGAL_API_H */
