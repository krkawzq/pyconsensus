/* hts_shim.c — thin accessors over htslib bcf1_t / bcf_hdr_t.
 *
 * Why this exists: bcf1_t has bitfields (n_allele:16) and a nested
 * bcf_dec_t (d.allele) whose exact layout is version-dependent.  Rather than
 * replicate that layout in Rust FFI (fragile), we expose a few trivial
 * accessors compiled against htslib's own headers.  This is the ONLY C glue
 * in the project; all consensus *algorithm* logic stays in Rust.
 *
 * Compiled by build.rs with: gcc -O2 -fPIC -I<libs/htslib> -c
 */
#include "htslib/vcf.h"
#include "htslib/hts.h"

long long shim_bcf_pos(const bcf1_t *r)      { return (long long)r->pos; }
long long shim_bcf_rlen(const bcf1_t *r)     { return (long long)r->rlen; }
int         shim_bcf_rid(const bcf1_t *r)    { return r->rid; }
int         shim_bcf_n_allele(const bcf1_t *r) { return r->n_allele; }

const char *shim_bcf_allele(const bcf1_t *r, int i) {
    if (i < 0 || i >= r->n_allele) return NULL;
    return r->d.allele[i];
}

int         shim_bcf_n_sample(const bcf1_t *r) { return r->n_sample; }

/* Wrappers for static-inline helpers so Rust can call them via FFI. */
int         shim_bcf_hdr_name2id(const bcf_hdr_t *h, const char *id) { return bcf_hdr_name2id(h, id); }
const char *shim_bcf_hdr_id2name(const bcf_hdr_t *h, int rid)        { return bcf_hdr_id2name(h, rid); }
int         shim_bcf_hdr_nsamples(const bcf_hdr_t *h)                { return bcf_hdr_nsamples(h); }
const char *shim_bcf_seqname(const bcf_hdr_t *h, const bcf1_t *r)    { return bcf_seqname(h, r); }
const char *shim_bcf_hdr_sample_name(const bcf_hdr_t *h, int i)      { return bcf_hdr_int2id(h, BCF_DT_SAMPLE, i); }
