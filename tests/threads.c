/*******************************************************************************
    threads: contexts on different threads don't interfere. Four threads,
    each with its own context and a 2 MiB stack (what callers such as Rust
    threads give the library), run metagenomic mode, single mode from a
    loaded model, and single mode trained on the input, several rounds each
    on the same context. Every round must give what the same run gave on the
    main thread (the first differing line is printed). Built with
    -fsanitize=thread, it also catches data races (the sanitizers CI job;
    docs/testing.md).

    usage: tests/threads   (from the repository root)
*******************************************************************************/

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prodigal.h"
#include "api_gff.h"
#include "fasta.h"

enum { META, MODEL, TRAIN, NKINDS };
static const char *kind_name[NKINDS] = { "meta", "model", "trained" };
static const int thread_kind[] = { META, MODEL, TRAIN, META };
#define NTHREADS ((int)(sizeof(thread_kind) / sizeof(thread_kind[0])))
#define ROUNDS 3
#define THREAD_STACK (2u << 20)

static fasta_t fa;
static void *model;
static size_t model_len;
static lines_t want[NKINDS];

/* One round of kind k on ctx: load or train the model, then find genes in
   every record, collecting the lines in out. 0, or -1 on an error. */
static int run(prodigal_ctx_t *ctx, int k, lines_t *out) {
    int rc = k == MODEL ? prodigal_load_training(ctx, model, model_len)
           : k == TRAIN ? api_train(ctx, &fa) : PRODIGAL_OK;
    if (rc == PRODIGAL_OK)
        rc = api_emit_records(ctx, &fa, 0, fa.n, k == META, lines_emit, out);
    if (rc != PRODIGAL_OK) {
        fprintf(stderr, "threads: %s: %s\n", kind_name[k], prodigal_last_error(ctx));
        return -1;
    }
    return 0;
}

/* Setup failures end the process (with any started threads). */
static void die(const char *what) {
    fprintf(stderr, "threads: %s\n", what);
    exit(1);
}

static prodigal_ctx_t *make(int k) {
    prodigal_config_t cfg;
    prodigal_config_init(&cfg);
    cfg.meta_mode = k == META;
    return prodigal_create(&cfg);
}

/* Thread body: ROUNDS rounds on one context; returns the number that
   differ from the main thread's run (as a pointer-sized integer). */
static void *worker(void *arg) {
    int k = thread_kind[(intptr_t)arg], r;
    intptr_t bad = 0;
    char what[64];
    prodigal_ctx_t *ctx = make(k);
    if (ctx == NULL) die("prodigal_create failed on a thread");
    for (r = 0; r < ROUNDS; r++) {
        lines_t got = {0};
        snprintf(what, sizeof(what), "thread %d (%s), round %d",
                 (int)(intptr_t)arg, kind_name[k], r + 1);
        bad += run(ctx, k, &got) != 0 || !lines_equal(&want[k], &got, what);
        lines_free(&got);
    }
    prodigal_destroy(ctx);
    return (void *)bad;
}

int main(void) {
    pthread_t tid[NTHREADS];
    pthread_attr_t attr;
    intptr_t i, bad = 0;
    void *ret;

    if (fasta_load("anthus_aco.fas", &fa) != 0 ||
        (model = load_file("testdata/ground_truth/ref_train.bin", &model_len)) == NULL)
        die("cannot read the fixtures (run from the repository root)");
    for (i = 0; i < NKINDS; i++) {
        prodigal_ctx_t *ctx = make((int)i);
        if (ctx == NULL) die("prodigal_create failed");
        if (run(ctx, (int)i, &want[i]) != 0) die("the serial run failed");
        prodigal_destroy(ctx);
    }

    if (pthread_attr_init(&attr) != 0 ||
        pthread_attr_setstacksize(&attr, THREAD_STACK) != 0)
        die("cannot set a 2 MiB thread stack");
    for (i = 0; i < NTHREADS; i++)
        if (pthread_create(&tid[i], &attr, worker, (void *)i) != 0)
            die("pthread_create failed");
    for (i = 0; i < NTHREADS; i++) {
        if (pthread_join(tid[i], &ret) != 0) die("pthread_join failed");
        bad += (intptr_t)ret;
    }
    pthread_attr_destroy(&attr);
    for (i = 0; i < NKINDS; i++) lines_free(&want[i]);
    fasta_free(&fa);
    free(model);
    printf("threads: %d threads x %d rounds %s\n", NTHREADS, ROUNDS,
           bad == 0 ? "match the serial runs" : "DIFFER");
    return bad != 0;
}
