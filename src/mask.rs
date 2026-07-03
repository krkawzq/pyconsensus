//! mask — `-m FILE` + `--mask-with CHAR|uc|lc`.
//!
//! Ports `mask_region` (consensus.c:1063) and the per-variant mask check
//! (consensus.c:590-600). Mask regions are loaded into an htslib regidx; at
//! apply time, variants overlapping a char-mode mask are skipped, and at the
//! end of the region the masked spans are overwritten per `--mask-with`.

use crate::htslib_ffi as ffi;
use std::ffi::CString;
use std::path::PathBuf;

/// `--mask-with` mode. Char-mode also causes variant skipping (MASK_SKIP).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaskWith {
    /// Replace masked bases with this char; skip overlapping variants.
    Char(u8),
    /// Uppercase masked bases; do NOT skip variants.
    Uc,
    /// Lowercase masked bases; do NOT skip variants.
    Lc,
}

impl Default for MaskWith {
    fn default() -> Self {
        MaskWith::Char(b'N')
    }
}

impl MaskWith {
    /// MASK_SKIP(mask) = (with != UC && with != LC) — char mode skips variants.
    pub fn skips_variants(&self) -> bool {
        !matches!(self, MaskWith::Uc | MaskWith::Lc)
    }
}

/// One mask file + its mode. Not Send (regidx/regitr carry state).
pub struct Mask {
    idx: *mut ffi::regidx_t,
    itr: *mut ffi::regitr_t,
    pub with: MaskWith,
    path: PathBuf,
}

impl Mask {
    /// Load a mask BED file with the given replacement mode.
    pub fn load(path: impl Into<PathBuf>, with: MaskWith) -> Result<Self, String> {
        let path = path.into();
        let c = CString::new(path.to_str().ok_or("non-UTF8 mask path")?)
            .map_err(|_| "non-NUL mask path".to_string())?;
        let idx = unsafe { ffi::regidx_init(c.as_ptr(), None, None, 0, std::ptr::null_mut()) };
        if idx.is_null() {
            return Err(format!("regidx_init failed for {}", path.display()));
        }
        let itr = unsafe { ffi::regitr_init(idx) };
        if itr.is_null() {
            unsafe { ffi::regidx_destroy(idx) };
            return Err("regitr_init failed".to_string());
        }
        Ok(Mask {
            idx,
            itr,
            with,
            path,
        })
    }

    pub fn path(&self) -> &std::path::Path {
        &self.path
    }

    /// Does any mask region overlap `[start, end]` (0-based, inclusive)?
    /// Used by apply_variant's MASK_SKIP check (consensus.c:590-600).
    pub fn overlaps(&self, chr: &str, start: i64, end: i64) -> bool {
        let c = match CString::new(chr) {
            Ok(c) => c,
            Err(_) => return false,
        };
        unsafe { ffi::regidx_overlap(self.idx, c.as_ptr(), start, end, self.itr) != 0 }
    }

    /// Apply the mask to `buf`, which spans `[ori_pos, ori_pos+len-1]`
    /// (0-based). Mirrors `mask_region` (consensus.c:1063): for each mask
    /// region overlapping the buffer, replace the overlapping span per `with`.
    pub fn apply_to_buf(&self, chr: &str, buf: &mut [u8], ori_pos: i64) {
        if buf.is_empty() {
            return;
        }
        let len = buf.len() as i64;
        // bcftools: start = fa_src_pos - len; end = fa_src_pos, where fa_src_pos
        // is the (0-based) position just past the buffer end. With our whole-
        // region buffer, fa_src_pos = ori_pos + len, so [start,end] = [ori_pos, ori_pos+len-1].
        let start = ori_pos;
        let end = ori_pos + len - 1;
        let c = match CString::new(chr) {
            Ok(c) => c,
            Err(_) => return,
        };
        unsafe {
            if ffi::regidx_overlap(self.idx, c.as_ptr(), start, end, self.itr) == 0 {
                return;
            }
            while ffi::regitr_overlap(self.itr) != 0 {
                let mbeg = (*self.itr).beg;
                let mend = (*self.itr).end;
                // map to buffer indices
                let mut idx_start = mbeg - start;
                let mut idx_end = mend - start;
                if idx_start < 0 {
                    idx_start = 0;
                }
                if idx_end >= len {
                    idx_end = len - 1;
                }
                if idx_end < idx_start {
                    continue;
                }
                let (s, e) = (idx_start as usize, idx_end as usize);
                match self.with {
                    MaskWith::Char(ch) => {
                        for b in &mut buf[s..=e] {
                            *b = ch;
                        }
                    }
                    MaskWith::Uc => {
                        for b in &mut buf[s..=e] {
                            *b = b.to_ascii_uppercase();
                        }
                    }
                    MaskWith::Lc => {
                        for b in &mut buf[s..=e] {
                            *b = b.to_ascii_lowercase();
                        }
                    }
                }
            }
        }
    }
}

impl Drop for Mask {
    fn drop(&mut self) {
        unsafe {
            if !self.itr.is_null() {
                ffi::regitr_destroy(self.itr);
            }
            if !self.idx.is_null() {
                ffi::regidx_destroy(self.idx);
            }
        }
    }
}

// regidx/regitr carry internal state; not shared across threads unsynchronised.
unsafe impl Send for Mask {}
