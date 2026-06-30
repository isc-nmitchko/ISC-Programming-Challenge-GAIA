/* gaiascan.c — hyper-optimized Gaia DR3 epoch-photometry kernel.
 *
 * Per file: decompress with libdeflate, then a single forward byte-scan that
 * parses ONLY source_id + the two flux arrays (bp_flux col 11, rp_flux col 16),
 * tracking per-band min/max running doubles. No DataFrame, no per-value alloc.
 * Files are processed in parallel with OpenMP; ctypes releases the GIL on the call.
 *
 * Build:
 *   gcc -O3 -fopenmp -shared -fPIC src/gaiascan.c -o src/gaiascan.so \
 *       /lib/x86_64-linux-gnu/libdeflate.so.0 -lm
 *
 * Result rows are returned to Python as a packed array of Rec; Python applies the
 * exact same percentage_change/filter/format as the polars path, so output is
 * byte-identical to the polars oracle.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* libdeflate has no -dev header in this image; declare the tiny ABI we use. */
struct libdeflate_decompressor;
extern struct libdeflate_decompressor *libdeflate_alloc_decompressor(void);
extern int libdeflate_gzip_decompress(struct libdeflate_decompressor *,
                                       const void *in, size_t in_nbytes,
                                       void *out, size_t out_nbytes_avail,
                                       size_t *actual_out_nbytes_ret);
extern void libdeflate_free_decompressor(struct libdeflate_decompressor *);

/* Scan a bracketed numeric array field [s,e) — track min/max of finite values.
 * NaN tokens and empty slots are skipped (strtod yields nan / no progress). */
static inline void scan_array(const char *s, const char *e, double *mn, double *mx) {
    double lo = INFINITY, hi = -INFINITY;
    char *q = (char *)s;
    while (q < e) {
        char c = *q;
        /* fast-skip the non-numeric framing: quotes, brackets, commas, spaces */
        if (c != '-' && c != '+' && c != '.' && (c < '0' || c > '9') &&
            c != 'n' && c != 'N' && c != 'i' && c != 'I') { q++; continue; }
        char *endp;
        double v = strtod(q, &endp);
        if (endp == q) { q++; continue; }
        q = endp;
        if (isfinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
    }
    *mn = (lo == INFINITY) ? NAN : lo;   /* no valid value -> NaN (Python maps to null) */
    *mx = (hi == -INFINITY) ? NAN : hi;
}

/* Read an entire file into a malloc'd buffer; returns size via *n (0 on error). */
static unsigned char *slurp(const char *path, size_t *n) {
    FILE *f = fopen(path, "rb");
    if (!f) { *n = 0; return NULL; }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *buf = (unsigned char *)malloc(sz);
    if (!buf) { fclose(f); *n = 0; return NULL; }
    size_t got = fread(buf, 1, sz, f);
    fclose(f);
    *n = got;
    return buf;
}

/* Output is struct-of-arrays so Python can bulk-copy each column with zero
 * per-element marshalling. ids[], and the four flux columns, filled at row idx. */
typedef struct {
    int64_t *ids;
    double *bp_min, *bp_max, *rp_min, *rp_max;
    long max_rows;
} Out;

/* Process one decompressed CSV buffer; emit one row per data line via atomic idx. */
static void scan_buffer(const char *buf, size_t len, Out *o, long *count) {
    const char *p = buf, *end = buf + len;
    while (p < end) {
        const char *nl = (const char *)memchr(p, '\n', end - p);
        if (!nl) nl = end;
        /* data rows begin with a digit (solution_id); skip '#' header and column header */
        if (*p >= '0' && *p <= '9') {
            int field = 0, inq = 0;
            const char *fs = p;
            int64_t sid = 0;
            double bmn = NAN, bmx = NAN, rmn = NAN, rmx = NAN;
            for (const char *c = p; c <= nl; c++) {
                if (c == nl || (*c == ',' && !inq)) {
                    if (field == 1) sid = strtoll(fs, NULL, 10);
                    else if (field == 11) scan_array(fs, c, &bmn, &bmx);
                    else if (field == 16) scan_array(fs, c, &rmn, &rmx);
                    field++;
                    fs = c + 1;
                    if (field > 16) break;   /* nothing useful past rp_flux */
                } else if (*c == '"') {
                    inq = !inq;
                }
            }
            long idx;
            #pragma omp atomic capture
            idx = (*count)++;
            if (idx < o->max_rows) {
                o->ids[idx] = sid;
                o->bp_min[idx] = bmn; o->bp_max[idx] = bmx;
                o->rp_min[idx] = rmn; o->rp_max[idx] = rmx;
            }
        }
        p = nl + 1;
    }
}

/* Public entry: decompress + scan every path in parallel. Returns row count.
 * Caller passes five preallocated column buffers (sized >= total rows). */
long gaia_scan(const char **paths, int npaths,
               int64_t *ids, double *bp_min, double *bp_max,
               double *rp_min, double *rp_max, long max_rows) {
    Out o = { ids, bp_min, bp_max, rp_min, rp_max, max_rows };
    long count = 0;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < npaths; i++) {
        size_t clen = 0;
        unsigned char *cbuf = slurp(paths[i], &clen);
        if (!cbuf || clen < 18) { free(cbuf); continue; }
        /* gzip trailer ISIZE = uncompressed size mod 2^32 (files are < 4 GB) */
        uint32_t isize = (uint32_t)cbuf[clen - 4] | ((uint32_t)cbuf[clen - 3] << 8) |
                         ((uint32_t)cbuf[clen - 2] << 16) | ((uint32_t)cbuf[clen - 1] << 24);
        unsigned char *ubuf = (unsigned char *)malloc(isize ? isize : 1);
        struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
        size_t actual = 0;
        int rc = libdeflate_gzip_decompress(d, cbuf, clen, ubuf, isize, &actual);
        libdeflate_free_decompressor(d);
        free(cbuf);
        if (rc == 0) scan_buffer((const char *)ubuf, actual, &o, &count);
        free(ubuf);
    }
    return count;
}
