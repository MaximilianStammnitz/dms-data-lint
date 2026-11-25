/*
 * dms_data_lint - integrity and format checks of paired-end FASTQ files, the raw-data gate of
 * nf-core/deepmutscan.
 *
 *   dms_data_lint [--max-base-qual FILE] R1 R2
 *
 * The checks and their messages are those of `fq lint` (St. Jude Rust Labs,
 * https://github.com/stjude-rust-labs/fq, MIT) in its default settings, version 0.12.0, so that
 * the two can be compared case by case; no code of fq is used. For each record pair, in this order
 * (R1 before R2 for each check):
 *   S003 name starts with '@'      S004 name, sequence, plus line and quality not empty
 *   S002 sequence in ACGTNacgtn    S001 plus line starts with '+'
 *   S005 sequence and quality of equal length
 *   S006 quality characters printable ASCII (33-126)
 *   P001 R1 and R2 names equal
 * and, once all of that passed, S007: no read name twice in R1. Names are compared after
 * truncation at their last '/' or ' '. A file is read as gzip when its name ends in ".gz"
 * (multi-member, CRC32 and length of every member checked), else as plain text. Lines end at '\n';
 * one '\r' before it is dropped. The first problem stops the check.
 *
 * Speed: ISA-L decompression in the process, one thread per file, vector instructions for the
 * character checks (GCC/Clang vector extensions: SSE2, AVX2 or NEON).
 *
 * S007: the 64-bit hash of every R1 name is kept (8 bytes a read) and sorted at the end; only if two
 * hashes are equal is R1 read a second time, to compare the names themselves. The result is exact;
 * without duplicates a second pass is needed with a probability of about n^2 / 2^65.
 *
 * Exit 0 and nothing on stdout when the pair passes. Otherwise exit 1 and one message on stderr:
 * "<file>:<line>:<column>: [CODE] Name: message" for a record, "Error: r1-src unexpectedly ended"
 * (or r2-src) for files with different numbers of records, "Error: I/O error in <file>: ..." for a
 * file that cannot be read or decompressed.
 * --max-base-qual FILE: for a pair that passes, the highest base quality (Phred+33) among the first
 * 10,000 records of each file, tab-separated.
 *
 * MIT License, see LICENSE.
 */
#define _GNU_SOURCE
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <isa-l.h>

#define DMS_DATA_LINT_VERSION "0.1.0"

#define CHUNK (4u << 20)
#define INBUF (1u << 20)
#define QDEPTH 4
#define QUAL_RECORDS 10000

/* ------------------------------------------------------------------ input */

typedef struct {
    int fd, gz, eof_in, done, pending_err;
    uint64_t nin;
    struct inflate_state st;
    uint8_t *in;
    const char *path;
    char err[256];
} src_t;

static int src_open(src_t *s, const char *path) {
    memset(s, 0, sizeof *s);
    s->path = path;
    s->fd = open(path, O_RDONLY);
    if (s->fd < 0) return -1;
    size_t n = strlen(path);
    s->gz = n > 3 && strcmp(path + n - 3, ".gz") == 0;
    if (s->gz) {
        s->in = malloc(INBUF);
        isal_inflate_init(&s->st);
        s->st.crc_flag = ISAL_GZIP;
    }
    return 0;
}

static int refill(src_t *s) {
    ssize_t r = read(s->fd, s->in, INBUF);
    if (r < 0) { snprintf(s->err, sizeof s->err, "%s", strerror(errno)); return -1; }
    if (r == 0) s->eof_in = 1;
    s->nin += (uint64_t)r;
    s->st.next_in = s->in;
    s->st.avail_in = (uint32_t)r;
    return 0;
}

static const char *isal_msg(int r) {
    switch (r) {
    case ISAL_INVALID_BLOCK: return "corrupt deflate stream (invalid block)";
    case ISAL_INVALID_SYMBOL: return "corrupt deflate stream (invalid symbol)";
    case ISAL_INVALID_LOOKBACK: return "corrupt deflate stream (invalid lookback distance)";
    case ISAL_INVALID_WRAPPER: return "invalid gzip header";
    case ISAL_UNSUPPORTED_METHOD: return "invalid gzip header (unsupported compression method)";
    case ISAL_INCORRECT_CHECKSUM: return "corrupt gzip stream does not have a matching checksum";
    default: return "corrupt deflate stream";
    }
}

/* Appends up to `room` decompressed bytes at `dst`. Returns bytes written; 0 at the end of the
 * input; -1 on a read error (message in s->err). */
static ssize_t src_read(src_t *s, uint8_t *dst, size_t room) {
    if (s->pending_err) return -1;   /* the bytes before a read error were handed out first */
    if (s->done) return 0;
    if (!s->gz) {
        ssize_t r = read(s->fd, dst, room);
        if (r < 0) { snprintf(s->err, sizeof s->err, "%s", strerror(errno)); return -1; }
        if (r == 0) s->done = 1;
        return r;
    }
    size_t got = 0;
    while (got < room) {
        if (s->st.avail_in == 0 && !s->eof_in && refill(s) < 0) return -1;
        if (s->st.block_state == ISAL_BLOCK_FINISH) {
            /* a member ended: the end of the file, or the next member's header */
            if (s->st.avail_in == 0) { if (s->eof_in) { s->done = 1; break; } continue; }
            uint8_t *ni = s->st.next_in; uint32_t na = s->st.avail_in;
            isal_inflate_reset(&s->st);
            s->st.crc_flag = ISAL_GZIP;
            s->st.next_in = ni; s->st.avail_in = na;
        }
        if (s->st.avail_in == 0 && s->eof_in) {   /* a member, or the first header, is incomplete */
            snprintf(s->err, sizeof s->err, "unexpected end of file");
            s->pending_err = 1;
            return got ? (ssize_t)got : -1;
        }
        s->st.next_out = dst + got;
        s->st.avail_out = (uint32_t)(room - got);
        int r = isal_inflate(&s->st);
        got = room - s->st.avail_out;
        if (r != ISAL_DECOMP_OK) {
            snprintf(s->err, sizeof s->err, "%s", isal_msg(r));
            s->pending_err = 1;
            return got ? (ssize_t)got : -1;
        }
    }
    return (ssize_t)got;
}

/* ------------------------------------------------------------------ records */

typedef struct { const uint8_t *p; uint32_t n; } line_t;

/* One record from [p, e): fq's reader takes four lines; missing lines are empty. Returns the
 * position after the record, or NULL if the record is not complete and more data may follow. */
static const uint8_t *take_record(const uint8_t *p, const uint8_t *e, int at_eof, line_t L[4]) {
    for (int k = 0; k < 4; k++) {
        if (p == e) {
            if (!at_eof) return NULL;
            L[k].p = p; L[k].n = 0; continue;
        }
        const uint8_t *nl = memchr(p, '\n', (size_t)(e - p));
        if (!nl) {
            if (!at_eof) return NULL;
            L[k].p = p; L[k].n = (uint32_t)(e - p); p = e; continue;
        }
        uint32_t n = (uint32_t)(nl - p);
        if (n && p[n - 1] == '\r') n--;
        L[k].p = p; L[k].n = n; p = nl + 1;
    }
    return p;
}

static uint32_t trunc_name(const line_t *nm) {
    for (uint32_t i = nm->n; i > 0; i--)
        if (nm->p[i - 1] == '/' || nm->p[i - 1] == ' ') return i - 1;
    return nm->n;
}

enum { V_S003, V_S004, V_S002, V_S001, V_S005, V_S006, NV };

typedef struct {
    uint64_t rec;        /* record index; UINT64_MAX: none */
    int mask;
    int s004_which;      /* 0 name, 1 sequence, 2 plus line, 3 quality */
    uint32_t s002_col, s006_col;
    uint8_t s002_chr, s006_chr;
    uint32_t seq_len, qual_len;
} fail_t;

/* x | 0x20 is one of acgtn exactly for ACGTNacgtn. GCC/Clang vector extensions: SSE2, AVX2 or
 * NEON, whatever the target has. */
typedef uint8_t v16 __attribute__((vector_size(16)));
static inline int seq_ok(const uint8_t *s, uint32_t n) {
    v16 bad = {0};
    uint32_t i = 0;
    for (; i + 16 <= n; i += 16) {
        v16 x; memcpy(&x, s + i, 16);
        x |= 0x20;
        bad |= ~((x == 'a') | (x == 'c') | (x == 'g') | (x == 't') | (x == 'n'));
    }
    uint8_t b = 0;
    for (int k = 0; k < 16; k++) b |= bad[k];
    for (; i < n; i++) { uint8_t x = s[i] | 0x20; b |= !(x == 'a' || x == 'c' || x == 'g' || x == 't' || x == 'n'); }
    return !b;
}
static inline int qual_ok(const uint8_t *q, uint32_t n) {
    uint8_t bad = 0;
    for (uint32_t i = 0; i < n; i++) bad |= (uint8_t)((uint8_t)(q[i] - 33) > 93);
    return !bad;
}
static inline uint8_t qual_max(const uint8_t *q, uint32_t n) {
    uint8_t m = 0;
    for (uint32_t i = 0; i < n; i++) m = q[i] > m ? q[i] : m;
    return m;
}
static int alpha_ok(uint8_t c) { uint8_t x = c | 0x20; return x == 'a' || x == 'c' || x == 'g' || x == 't' || x == 'n'; }

/* fq's single-read validators on one record; bits in validator order */
static int check_record(const line_t L[4], uint32_t tl, fail_t *f) {
    int m = 0;
    if (tl == 0 || L[0].p[0] != '@') m |= 1 << V_S003;
    if (tl == 0) { m |= 1 << V_S004; f->s004_which = 0; }
    else if (L[1].n == 0) { m |= 1 << V_S004; f->s004_which = 1; }
    else if (L[2].n == 0) { m |= 1 << V_S004; f->s004_which = 2; }
    else if (L[3].n == 0) { m |= 1 << V_S004; f->s004_which = 3; }
    if (!seq_ok(L[1].p, L[1].n)) {
        m |= 1 << V_S002;
        for (uint32_t i = 0; i < L[1].n; i++) if (!alpha_ok(L[1].p[i])) { f->s002_col = i + 1; f->s002_chr = L[1].p[i]; break; }
    }
    if (L[2].n == 0 || L[2].p[0] != '+') m |= 1 << V_S001;
    if (L[1].n != L[3].n) { m |= 1 << V_S005; f->seq_len = L[1].n; f->qual_len = L[3].n; }
    if (!qual_ok(L[3].p, L[3].n)) {
        m |= 1 << V_S006;
        for (uint32_t i = 0; i < L[3].n; i++) if ((uint8_t)(L[3].p[i] - 33) > 93) { f->s006_col = i + 1; f->s006_chr = L[3].p[i]; break; }
    }
    return m;
}

static inline uint64_t hash_name(const uint8_t *p, uint32_t n) {
    uint64_t h = 0x9E3779B97F4A7C15ull ^ n;
    while (n >= 8) {
        uint64_t w; memcpy(&w, p, 8);
        h = (h ^ w) * 0xff51afd7ed558ccdull; h ^= h >> 32;
        p += 8; n -= 8;
    }
    uint64_t w = 0; memcpy(&w, p, n);
    h = (h ^ w) * 0xc4ceb9fe1a85ec53ull;
    h ^= h >> 33; h *= 0xff51afd7ed558ccdull; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ull; h ^= h >> 33;
    return h;
}

/* ------------------------------------------------------------------ chunks and queues */

typedef struct {
    uint8_t *buf; size_t cap, len;
    uint32_t *noff, *nlen; size_t nrec, rcap;
    int last;            /* final chunk of the stream */
    int ioerr;           /* the stream stopped at a read error after this chunk */
} chunk_t;

typedef struct {
    chunk_t *q[QDEPTH]; int head, n;
    chunk_t *pool[QDEPTH + 2]; int npool;
    pthread_mutex_t mu; pthread_cond_t cv;
} queue_t;

static void q_push(queue_t *q, chunk_t *c) {
    pthread_mutex_lock(&q->mu);
    while (q->n == QDEPTH) pthread_cond_wait(&q->cv, &q->mu);
    q->q[(q->head + q->n) % QDEPTH] = c; q->n++;
    pthread_cond_broadcast(&q->cv);
    pthread_mutex_unlock(&q->mu);
}
static chunk_t *q_pop(queue_t *q) {
    pthread_mutex_lock(&q->mu);
    while (q->n == 0) pthread_cond_wait(&q->cv, &q->mu);
    chunk_t *c = q->q[q->head]; q->head = (q->head + 1) % QDEPTH; q->n--;
    pthread_cond_broadcast(&q->cv);
    pthread_mutex_unlock(&q->mu);
    return c;
}
static void q_recycle(queue_t *q, chunk_t *c) {
    pthread_mutex_lock(&q->mu);
    if (q->npool < QDEPTH + 2) q->pool[q->npool++] = c; else { free(c->buf); free(c->noff); free(c->nlen); free(c); }
    pthread_mutex_unlock(&q->mu);
}
static chunk_t *q_fresh(queue_t *q) {
    chunk_t *c = NULL;
    pthread_mutex_lock(&q->mu);
    if (q->npool) c = q->pool[--q->npool];
    pthread_mutex_unlock(&q->mu);
    if (!c) {
        c = calloc(1, sizeof *c);
        c->cap = CHUNK + (1u << 16); c->buf = malloc(c->cap);
        c->rcap = 1u << 15; c->noff = malloc(c->rcap * 4); c->nlen = malloc(c->rcap * 4);
    }
    c->len = 0; c->nrec = 0; c->last = 0; c->ioerr = 0;
    return c;
}

/* ------------------------------------------------------------------ per-file worker */

typedef struct {
    src_t src;
    queue_t q;
    fail_t fail;
    uint64_t nrec;
    uint8_t qmax;
    int hash_names;      /* R1: keep every name's hash for S007 */
    /* R1 name hashes in 256 buckets by their top byte, each a list of fixed blocks, so that memory
     * stays at 8 bytes a read (no doubling) */
    struct hblk { struct hblk *next; size_t n; uint64_t h[8190]; } *bk[256];
    size_t bn[256];
} worker_t;

static void keep_hash(worker_t *w, uint64_t h) {
    unsigned b = (unsigned)(h >> 56);
    struct hblk *k = w->bk[b];
    if (!k || k->n == 8190) { struct hblk *x = malloc(sizeof *x); x->next = k; x->n = 0; w->bk[b] = k = x; }
    k->h[k->n++] = h; w->bn[b]++;
}

static void *worker(void *arg) {
    worker_t *w = arg;
    uint8_t *carry = NULL; size_t ncarry = 0, ccap = 0;
    for (;;) {
        chunk_t *c = q_fresh(&w->q);
        size_t target = ncarry < CHUNK / 2 ? CHUNK : ncarry * 2;   /* a record longer than a chunk doubles it */
        if (c->cap < target + (1u << 16)) { c->cap = target + (1u << 16); c->buf = realloc(c->buf, c->cap); }
        memcpy(c->buf, carry, ncarry); c->len = ncarry; ncarry = 0;
        int at_eof = 0, err = 0;
        while (c->len < target) {
            ssize_t r = src_read(&w->src, c->buf + c->len, c->cap - c->len);
            if (r < 0) { err = 1; break; }
            if (r == 0) { at_eof = 1; break; }
            c->len += (size_t)r;
        }
        const uint8_t *p = c->buf, *e = c->buf + c->len;
        for (;;) {
            line_t L[4];
            if (p == e && (at_eof || err)) break;
            const uint8_t *nx = take_record(p, e, at_eof, L);
            if (!nx) break;
            uint32_t tl = trunc_name(&L[0]);
            if (w->fail.rec == UINT64_MAX) {
                fail_t f = {0};
                int m = check_record(L, tl, &f);
                if (m) { f.rec = w->nrec; f.mask = m; w->fail = f; }
            }
            if (w->nrec < QUAL_RECORDS) { uint8_t m = qual_max(L[3].p, L[3].n); if (m > w->qmax) w->qmax = m; }
            if (w->hash_names) keep_hash(w, hash_name(L[0].p, tl));
            if (c->nrec == c->rcap) { c->rcap *= 2; c->noff = realloc(c->noff, c->rcap * 4); c->nlen = realloc(c->nlen, c->rcap * 4); }
            c->noff[c->nrec] = (uint32_t)(L[0].p - c->buf); c->nlen[c->nrec] = tl; c->nrec++;
            w->nrec++;
            p = nx;
        }
        if (!at_eof && !err) {
            ncarry = (size_t)(e - p);
            if (ncarry > ccap) { ccap = ncarry * 2; carry = realloc(carry, ccap); }
            memcpy(carry, p, ncarry);
        }
        c->last = at_eof || err; c->ioerr = err;
        q_push(&w->q, c);
        if (c->last) break;
    }
    free(carry);
    return NULL;
}

/* ------------------------------------------------------------------ reporting */

static void put_char(FILE *o, uint8_t c) {      /* Rust's char::from(u8): Latin-1 code point */
    if (c < 0x80) fputc(c, o); else { fputc(0xC0 | (c >> 6), o); fputc(0x80 | (c & 0x3F), o); }
}
static void put_lossy(FILE *o, const uint8_t *p, uint32_t n) { fwrite(p, 1, n, o); }

static const char *VCODE[NV] = {"S003", "S004", "S002", "S001", "S005", "S006"};
static const char *VNAME[NV] = {"NameValidator", "CompleteValidator", "AlphabetValidator", "PlusLineValidator", "ConsistentSeqQualValidator", "QualityStringValidator"};

static void report(const char *path, uint64_t rec, int v, const fail_t *f) {
    static const int line_of_s004[4] = {0, 1, 2, 3};
    int lt; uint32_t col = 1;
    switch (v) {
    case V_S003: lt = 0; break;
    case V_S004: lt = line_of_s004[f->s004_which]; break;
    case V_S002: lt = 1; col = f->s002_col; break;
    case V_S001: lt = 2; break;
    case V_S005: lt = 1; break;
    default: lt = 3; col = f->s006_col; break;
    }
    fprintf(stderr, "%s:%llu:%u: [%s] %s: ", path, (unsigned long long)(rec * 4 + lt + 1), col, VCODE[v], VNAME[v]);
    static const char *s004[4] = {"empty name", "empty sequence", "empty plus line", "empty quality"};
    switch (v) {
    case V_S003: fputs("missing @ prefix", stderr); break;
    case V_S004: fputs(s004[f->s004_which], stderr); break;
    case V_S002: fputs("invalid sequence character: '", stderr); put_char(stderr, f->s002_chr); fputc('\'', stderr); break;
    case V_S001: fputs("missing + prefix", stderr); break;
    case V_S005: fprintf(stderr, "sequence-quality scores lengths mismatch: expected %u, got %u", f->seq_len, f->qual_len); break;
    default: fputs("invalid quality score character: '", stderr); put_char(stderr, f->s006_chr); fputc('\'', stderr); break;
    }
    fputc('\n', stderr);
}

static int cmp_u64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }

static void radix56(uint64_t *a, size_t n, uint64_t *tmp) {   /* sorts on bits 0-55 */
    static size_t cnt[1 << 14];
    int bits = n > 65536 ? 14 : 8, passes = (56 + bits - 1) / bits, nb = 1 << bits;
    uint64_t mask = (uint64_t)nb - 1;
    for (int pass = 0; pass < passes; pass++) {
        int sh = pass * bits;
        memset(cnt, 0, (size_t)nb * sizeof cnt[0]);
        for (size_t i = 0; i < n; i++) cnt[(a[i] >> sh) & mask]++;
        size_t s = 0; for (int i = 0; i < nb; i++) { size_t t = cnt[i]; cnt[i] = s; s += t; }
        for (size_t i = 0; i < n; i++) tmp[cnt[(a[i] >> sh) & mask]++] = a[i];
        memcpy(a, tmp, n * 8);
    }
}

/* S007, second pass: R1 again, exact names for the hashes seen twice */
typedef struct { uint8_t *name; uint32_t n; int seen; } slot_t;

static int s007_confirm(const char *path, uint64_t *sus, size_t nsus) {
    size_t cap = 16; while (cap < nsus * 4) cap <<= 1;
    slot_t *tab = calloc(cap, sizeof *tab);
    src_t s; if (src_open(&s, path) < 0) { fprintf(stderr, "Error: I/O error in %s: %s\n", path, strerror(errno)); return 2; }
    size_t bcap = CHUNK + (1u << 16), len = 0; uint8_t *buf = malloc(bcap);
    uint64_t rec = 0; int at_eof = 0, ret = 0;
    for (;;) {
        while (!at_eof && len < CHUNK) {
            ssize_t r = src_read(&s, buf + len, bcap - len);
            if (r < 0) { fprintf(stderr, "Error: I/O error in %s: %s\n", path, s.err); ret = 2; goto out; }
            if (r == 0) at_eof = 1; else len += (size_t)r;
        }
        const uint8_t *p = buf, *e = buf + len;
        for (;;) {
            line_t L[4];
            if (p == e && at_eof) break;
            const uint8_t *nx = take_record(p, e, at_eof, L);
            if (!nx) break;
            uint32_t tl = trunc_name(&L[0]);
            uint64_t h = hash_name(L[0].p, tl);
            if (bsearch(&h, sus, nsus, 8, cmp_u64)) {
                size_t i = (size_t)(h & (cap - 1));
                for (;; i = (i + 1) & (cap - 1)) {
                    if (!tab[i].name) { tab[i].name = malloc(tl ? tl : 1); memcpy(tab[i].name, L[0].p, tl); tab[i].n = tl; break; }
                    if (tab[i].n == tl && memcmp(tab[i].name, L[0].p, tl) == 0) {
                        fprintf(stderr, "%s:%llu:1: [S007] DuplicateNameValidator: duplicate name: '", path, (unsigned long long)(rec * 4 + 1));
                        put_lossy(stderr, L[0].p, tl); fputs("'\n", stderr);
                        ret = 1; goto out;
                    }
                }
            }
            rec++; p = nx;
        }
        if (at_eof) break;
        len = (size_t)(e - p); memmove(buf, p, len);
        if (len > bcap / 2) { bcap *= 2; buf = realloc(buf, bcap); }
    }
out:
    free(buf);
    return ret;
}

/* ------------------------------------------------------------------ main */

static worker_t W[2];

static void usage(FILE *o) {
    fputs("dms_data_lint " DMS_DATA_LINT_VERSION " - integrity and format checks of paired-end FASTQ files\n\n"
          "usage: dms_data_lint [--max-base-qual FILE] R1 R2\n\n"
          "  R1, R2               the two files of a read pair, plain or gzip-compressed (.gz)\n"
          "  --max-base-qual FILE write the highest base quality of the first 10,000 records\n"
          "                       of R1 and R2 (Phred, tab-separated) for a pair that passes\n"
          "  --version            print the version of dms_data_lint and ISA-L\n"
          "  --help               print this help\n\n"
          "Exit status 0: the pair passes; 1: it does not (the reason on stderr); 2: usage error.\n", o);
}

int main(int argc, char **argv) {
    const char *qout = NULL; int a = 1;
    if (argc == 2 && strcmp(argv[1], "--version") == 0) {
#ifdef ISAL_MAJOR_VERSION
        printf("dms_data_lint %s (isa-l %d.%d.%d)\n", DMS_DATA_LINT_VERSION, ISAL_MAJOR_VERSION, ISAL_MINOR_VERSION, ISAL_PATCH_VERSION);
#else
        printf("dms_data_lint %s (isa-l)\n", DMS_DATA_LINT_VERSION);
#endif
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) { usage(stdout); return 0; }
    if (a < argc && strcmp(argv[a], "--max-base-qual") == 0 && a + 1 < argc) { qout = argv[a + 1]; a += 2; }
    if (argc - a != 2) { usage(stderr); return 2; }
    const char *path[2] = {argv[a], argv[a + 1]};
    for (int i = 0; i < 2; i++) {
        if (src_open(&W[i].src, path[i]) < 0) { fprintf(stderr, "Error: I/O error in %s: %s\n", path[i], strerror(errno)); return 1; }
        W[i].fail.rec = UINT64_MAX; W[i].hash_names = (i == 0);
        pthread_mutex_init(&W[i].q.mu, NULL); pthread_cond_init(&W[i].q.cv, NULL);
    }
    pthread_t th[2];
    for (int i = 0; i < 2; i++) pthread_create(&th[i], NULL, worker, &W[i]);

    chunk_t *c[2] = {NULL, NULL}; size_t ix[2] = {0, 0}; int ended[2] = {0, 0};
    uint64_t k = 0; int status = 0;
    for (;; k++) {
        int have[2];
        for (int i = 0; i < 2; i++) {
            while (!ended[i] && (!c[i] || ix[i] == c[i]->nrec)) {
                if (c[i] && c[i]->last) { ended[i] = c[i]->ioerr ? 2 : 1; break; }
                if (c[i]) q_recycle(&W[i].q, c[i]);
                c[i] = q_pop(&W[i].q); ix[i] = 0;
            }
            have[i] = !ended[i];
        }
        /* fq reads R1's record, then R2's, then compares the ends of the files */
        int stop = 0;
        for (int i = 0; i < 2 && !stop; i++)
            if (!have[i] && ended[i] == 2) { fprintf(stderr, "Error: I/O error in %s: %s\n", path[i], W[i].src.err); status = 1; stop = 1; }
        if (stop) break;
        if (!have[0] && !have[1]) break;
        if (!have[0] || !have[1]) { fprintf(stderr, "Error: %s unexpectedly ended\n", have[0] ? "r2-src" : "r1-src"); status = 1; break; }
        const uint8_t *n0 = c[0]->buf + c[0]->noff[ix[0]], *n1 = c[1]->buf + c[1]->noff[ix[1]];
        uint32_t l0 = c[0]->nlen[ix[0]], l1 = c[1]->nlen[ix[1]];
        int m0 = W[0].fail.rec == k ? W[0].fail.mask : 0, m1 = W[1].fail.rec == k ? W[1].fail.mask : 0;
        int p001 = l0 != l1 || memcmp(n0, n1, l0) != 0;
        if (m0 | m1 | p001) {
            for (int v = 0; v < NV; v++) {
                if (m0 & (1 << v)) { report(path[0], k, v, &W[0].fail); status = 1; break; }
                if (m1 & (1 << v)) { report(path[1], k, v, &W[1].fail); status = 1; break; }
            }
            if (!status) {
                fprintf(stderr, "%s:%llu:1: [P001] NamesValidator: names mismatch: expected '", path[0], (unsigned long long)(k * 4 + 1));
                put_lossy(stderr, n1, l1); fputs("', got '", stderr); put_lossy(stderr, n0, l0); fputs("'\n", stderr);
                status = 1;
            }
            break;
        }
        ix[0]++; ix[1]++;
    }
    if (status) _exit(1);   /* fq stops at the first error; the readers need not finish */
    for (int i = 0; i < 2; i++) pthread_join(th[i], NULL);

    /* S007: R1 names whose hashes are equal */
    size_t nmax = 0; for (int b = 0; b < 256; b++) if (W[0].bn[b] > nmax) nmax = W[0].bn[b];
    uint64_t *v = malloc((nmax ? nmax : 1) * 8), *tmp = malloc((nmax ? nmax : 1) * 8), *sus = NULL; size_t nsus = 0, scap = 0;
    for (int b = 0; b < 256; b++) {
        size_t n = 0;
        for (struct hblk *k = W[0].bk[b], *nx; k; k = nx) { memcpy(v + n, k->h, k->n * 8); n += k->n; nx = k->next; free(k); }
        radix56(v, n, tmp);
        for (size_t i = 1; i < n; i++) if (v[i] == v[i - 1] && (nsus == 0 || sus[nsus - 1] != v[i])) {
            if (nsus == scap) { scap = scap ? scap * 2 : 64; sus = realloc(sus, scap * 8); }
            sus[nsus++] = v[i];
        }
    }
    free(v); free(tmp);
    if (nsus) {
        qsort(sus, nsus, 8, cmp_u64);
        if (s007_confirm(path[0], sus, nsus)) return 1;
    }
    if (qout) {
        FILE *o = fopen(qout, "w");
        if (!o) { fprintf(stderr, "Error: could not create %s: %s\n", qout, strerror(errno)); return 1; }
        fprintf(o, "%d\t%d\n", W[0].qmax > 33 ? W[0].qmax - 33 : 0, W[1].qmax > 33 ? W[1].qmax - 33 : 0);
        fclose(o);
    }
    fprintf(stderr, "dms_data_lint: %llu read pairs, PASS\n", (unsigned long long)k);
    return 0;
}
