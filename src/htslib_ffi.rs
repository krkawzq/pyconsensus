//! Minimal, non-intrusive `extern "C"` bindings to htslib.
//!
//! Design (docs/design.md §5.5): we hand-declare only the htslib C functions
//! this project uses, instead of pulling in `rust-htslib`/`noodles`. htslib
//! macros / inline helpers are replicated in Rust (the GT decoding below, the
//! iupac tables in `iupac.rs`) rather than declared as FFI.
//!
//! Struct-field accessors that would require replicating version-dependent
//! layout (bcf1_t bitfields, nested bcf_dec_t) live in the tiny C shim
//! `src/hts_shim.c` and are declared here as `shim_*`.

use std::os::raw::{c_char, c_int, c_void};

// ===========================================================================
// Opaque htslib types
// ===========================================================================

/// Opaque htslib faidx handle.
#[repr(C)]
pub struct faidx_t {
    _private: [u8; 0],
}
#[repr(C)]
pub struct htsFile {
    _private: [u8; 0],
}
#[repr(C)]
pub struct bcf_hdr_t {
    _private: [u8; 0],
}
#[repr(C)]
pub struct bcf1_t {
    _private: [u8; 0],
}

/// Opaque htslib regidx handle (for `-m` mask).
#[repr(C)]
pub struct regidx_t {
    _private: [u8; 0],
}

/// htslib regitr_t — public struct (htslib/regidx.h:83). We only read beg/end.
#[repr(C)]
pub struct regitr_t {
    pub beg: i64,
    pub end: i64,
    pub payload: *mut c_void,
    pub seq: *mut c_char,
    pub itr: *mut c_void,
}

// ===========================================================================
// Constants (replicated from htslib headers)
// ===========================================================================

/// `enum fai_load_options { FAI_CREATE = 0x01 };`
pub const FAI_CREATE: c_int = 0x01;

/// bcf_unpack `which` flags
pub const BCF_UN_STR: c_int = 1;
pub const BCF_UN_FMT: c_int = 8;
pub const BCF_UN_ALL: c_int = BCF_UN_STR | 2 | 4 | BCF_UN_FMT;

/// bcf_hdr_id2int `type`
pub const BCF_DT_ID: c_int = 0;
pub const BCF_DT_CTG: c_int = 1;

/// bcf_get_format_values `type`
pub const BCF_HT_INT: c_int = 1;

/// Variant type bitmask (bcf_get_variant_types)
pub const VCF_REF: c_int = 0;
pub const VCF_SNP: c_int = 1 << 0;
pub const VCF_MNP: c_int = 1 << 1;
pub const VCF_INDEL: c_int = 1 << 2;
pub const VCF_OTHER: c_int = 1 << 3;
pub const VCF_BND: c_int = 1 << 4;
pub const VCF_OVERLAP: c_int = 1 << 5;

/// BCF int32 sentinels (for the flat GT array returned by bcf_get_format_values).
pub const BCF_INT32_VECTOR_END: i32 = i32::MIN + 1; // -2147483647
pub const BCF_INT32_MISSING: i32 = i32::MIN; // -2147483648

/// htslib log levels (`htslib/hts_log.h`).
pub const HTS_LOG_OFF: c_int = 0;
pub const HTS_LOG_ERROR: c_int = 1;
pub const HTS_LOG_WARNING: c_int = 3;
pub const HTS_LOG_INFO: c_int = 4;
pub const HTS_LOG_DEBUG: c_int = 5;
pub const HTS_LOG_TRACE: c_int = 6;

// ===========================================================================
// GT decoding (replicated from htslib/vcf.h bcf_gt_* macros)
// ===========================================================================
//
// Each allele in the GT array is an encoded int32:
//   bcf_gt_phased(idx)   = ((idx)+1) << 1 | 1
//   bcf_gt_unphased(idx) = ((idx)+1) << 1
//   bcf_gt_missing       = 0            (low bit may still carry phase info)
// so:
//   missing  <=> (raw >> 1) == 0   (i.e. raw in {0, 1})
//   allele   =  (raw >> 1) - 1     (-1 when missing; 0=REF, 1..=ALT idx)
//   phased   =  (raw & 1) != 0

#[inline]
pub fn gt_is_missing(raw: i32) -> bool {
    (raw >> 1) == 0
}

#[inline]
pub fn gt_allele(raw: i32) -> i32 {
    (raw >> 1) - 1
}

#[inline]
pub fn gt_is_phased(raw: i32) -> bool {
    (raw & 1) != 0
}

// ===========================================================================
// faidx
// ===========================================================================

extern "C" {
    pub fn fai_load(fn_: *const c_char) -> *mut faidx_t;
    pub fn fai_load3(
        fn_: *const c_char,
        fnfai: *const c_char,
        fngzi: *const c_char,
        flags: c_int,
    ) -> *mut faidx_t;
    pub fn fai_destroy(fai: *mut faidx_t);

    /// Fetch `[p_beg_i, p_end_i]` (0-based, inclusive). malloc'd; caller frees.
    pub fn faidx_fetch_seq64(
        fai: *const faidx_t,
        c_name: *const c_char,
        p_beg_i: i64,
        p_end_i: i64,
        len: *mut i64,
    ) -> *mut c_char;

    pub fn faidx_nseq(fai: *const faidx_t) -> c_int;
    pub fn faidx_has_seq(fai: *const faidx_t, seq: *const c_char) -> c_int;
    pub fn faidx_seq_len64(fai: *const faidx_t, seq: *const c_char) -> i64;

    pub fn free(ptr: *mut c_void);
}

// ===========================================================================
// hts / vcf / bcf (real, non-inline functions)
// ===========================================================================

extern "C" {
    pub fn hts_set_log_level(level: c_int);
    pub fn hts_get_log_level() -> c_int;

    pub fn hts_open(fn_: *const c_char, mode: *const c_char) -> *mut htsFile;
    pub fn hts_close(fp: *mut htsFile) -> c_int;

    pub fn bcf_hdr_read(fp: *mut htsFile) -> *mut bcf_hdr_t;
    pub fn bcf_hdr_destroy(h: *mut bcf_hdr_t);
    pub fn bcf_hdr_append(h: *mut bcf_hdr_t, line: *const c_char) -> c_int;
    pub fn bcf_hdr_sync(h: *mut bcf_hdr_t) -> c_int;
    /// `type` = BCF_DT_CTG for contigs. Returns -1 if not found.
    pub fn bcf_hdr_id2int(h: *const bcf_hdr_t, typ: c_int, id: *const c_char) -> c_int;
    /// Fills *nseqs; returns array of seq names (owned by header, do not free).
    pub fn bcf_hdr_seqnames(h: *const bcf_hdr_t, nseqs: *mut c_int) -> *mut *const c_char;

    pub fn bcf_init() -> *mut bcf1_t;
    /// Returns >=0 on success, -1 on EOF, <-1 on error.
    pub fn bcf_read(fp: *mut htsFile, h: *const bcf_hdr_t, v: *mut bcf1_t) -> c_int;
    pub fn bcf_destroy(v: *mut bcf1_t);
    pub fn bcf_unpack(b: *mut bcf1_t, which: c_int) -> c_int;

    /// Underlying impl of the `bcf_get_genotypes` macro.
    /// `tag` = "GT", `type` = BCF_HT_INT. Returns ploidy (values per sample)
    /// on success; negative on error. `*dst` is (re)allocated, caller frees.
    pub fn bcf_get_format_values(
        hdr: *const bcf_hdr_t,
        line: *mut bcf1_t,
        tag: *const c_char,
        dst: *mut *mut c_void,
        ndst: *mut c_int,
        typ: c_int,
    ) -> c_int;

    /// Bitmask of variant types across all ALTs.
    pub fn bcf_get_variant_types(rec: *mut bcf1_t) -> c_int;
    /// Variant type of the `ith_allele`-th ALT.
    pub fn bcf_get_variant_type(rec: *mut bcf1_t, ith_allele: c_int) -> c_int;
}

// ===========================================================================
// regidx (mask regions, `-m`)
// ===========================================================================

extern "C" {
    /// Load a region file (BED/tab) into an index. parsef=NULL → auto BED/tab.
    pub fn regidx_init(
        fname: *const c_char,
        parsef: Option<unsafe extern "C" fn() -> c_int>,
        freef: Option<unsafe extern "C" fn()>,
        payload_size: usize,
        usr: *mut c_void,
    ) -> *mut regidx_t;
    pub fn regidx_destroy(idx: *mut regidx_t);
    /// Query overlap with `[beg, end]` (0-based, inclusive). Returns 0 if no
    /// overlap, >0 if overlap (itr then iterates matches via regitr_overlap).
    pub fn regidx_overlap(
        idx: *mut regidx_t,
        chr: *const c_char,
        beg: i64,
        end: i64,
        itr: *mut regitr_t,
    ) -> c_int;
    pub fn regitr_init(idx: *mut regidx_t) -> *mut regitr_t;
    pub fn regitr_destroy(itr: *mut regitr_t);
    /// Advance itr to the next overlapping region. Returns 0 when exhausted.
    pub fn regitr_overlap(itr: *mut regitr_t) -> c_int;
}

// ===========================================================================
// Shim accessors (src/hts_shim.c) — avoid fragile struct layout in Rust
// ===========================================================================

extern "C" {
    pub fn shim_bcf_pos(r: *const bcf1_t) -> i64;
    pub fn shim_bcf_rlen(r: *const bcf1_t) -> i64;
    pub fn shim_bcf_rid(r: *const bcf1_t) -> c_int;
    pub fn shim_bcf_n_allele(r: *const bcf1_t) -> c_int;
    pub fn shim_bcf_n_sample(r: *const bcf1_t) -> c_int;
    /// Returns the i-th allele string (REF=0, ALTs=1..); NULL if out of range.
    pub fn shim_bcf_allele(r: *const bcf1_t, i: c_int) -> *const c_char;

    pub fn shim_bcf_hdr_name2id(h: *const bcf_hdr_t, id: *const c_char) -> c_int;
    pub fn shim_bcf_hdr_id2name(h: *const bcf_hdr_t, rid: c_int) -> *const c_char;
    pub fn shim_bcf_hdr_nsamples(h: *const bcf_hdr_t) -> c_int;
    pub fn shim_bcf_seqname(h: *const bcf_hdr_t, r: *const bcf1_t) -> *const c_char;
    pub fn shim_bcf_hdr_sample_name(h: *const bcf_hdr_t, i: c_int) -> *const c_char;
}

// ===========================================================================
// Thin safe wrapper: faidx handle (M0)
// ===========================================================================

/// An owned faidx handle. Drops via `fai_destroy`.
pub struct FaidxHandle {
    fai: *mut faidx_t,
}

impl FaidxHandle {
    /// Load a FASTA (building `.fai` if absent). Returns `None` on failure.
    pub fn load(path: &str) -> Option<Self> {
        let c = std::ffi::CString::new(path).ok()?;
        let fai = unsafe { fai_load(c.as_ptr()) };
        if fai.is_null() {
            None
        } else {
            Some(FaidxHandle { fai })
        }
    }

    pub fn nseq(&self) -> i32 {
        unsafe { faidx_nseq(self.fai) as i32 }
    }

    pub fn has_seq(&self, name: &str) -> bool {
        let c = match std::ffi::CString::new(name) {
            Ok(c) => c,
            Err(_) => return false,
        };
        unsafe { faidx_has_seq(self.fai, c.as_ptr()) != 0 }
    }

    pub fn seq_len(&self, name: &str) -> Option<i64> {
        let c = std::ffi::CString::new(name).ok()?;
        let n = unsafe { faidx_seq_len64(self.fai, c.as_ptr()) };
        if n < 0 {
            None
        } else {
            Some(n)
        }
    }

    /// Fetch `[beg, end]` (0-based, inclusive) as owned bytes.
    pub fn fetch(&self, name: &str, beg: i64, end: i64) -> Result<Vec<u8>, String> {
        let c = std::ffi::CString::new(name).map_err(|_| "invalid seq name".to_string())?;
        let mut len: i64 = 0;
        let p = unsafe { faidx_fetch_seq64(self.fai, c.as_ptr(), beg, end, &mut len) };
        if p.is_null() {
            return Err(format!("faidx_fetch_seq64 returned null (len={})", len));
        }
        if len < 0 {
            unsafe { free(p as *mut c_void) };
            return Err(format!("faidx_fetch_seq64 error len={}", len));
        }
        let slice = unsafe { std::slice::from_raw_parts(p as *const u8, len as usize) };
        let out = slice.to_vec();
        unsafe { free(p as *mut c_void) };
        Ok(out)
    }
}

impl Drop for FaidxHandle {
    fn drop(&mut self) {
        unsafe { fai_destroy(self.fai) }
    }
}

// faidx_t may carry mutable internal cursor state. M1's RefIndex hands out
// per-thread handles instead of sharing one across threads.
unsafe impl Send for FaidxHandle {}
