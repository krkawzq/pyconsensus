//! RefIndex — faidx wrapper, per-region fetch.
//!
//! (docs/design.md §5.3) The reference FASTA is "raw material" too, but its
//! preprocessing (building the `.fai` index) is cheap, so unlike VCF we do NOT
//! read the whole genome into memory. We keep a faidx handle and fetch each
//! requested region on demand into a reusable buffer.
//!
//! Coordinate convention (docs/implementation_plan.md §3):
//!   * `faidx_fetch_seq64` takes a **0-based, inclusive, closed** interval.
//!   * The Python API / script pass **1-based** regions; convert at the boundary.

use crate::htslib_ffi::FaidxHandle;
use crate::logging::ensure_default_htslib_log_level;
use std::path::PathBuf;
use std::sync::Mutex;

/// Owns the reference FASTA path plus a faidx handle.
///
/// `faidx_t` carries mutable cursor state, so the handle is guarded by a Mutex
/// to make `RefIndex: Send + Sync`. Contention is bounded: each fetch holds the
/// lock only for the duration of `faidx_fetch_seq64` (microseconds).
pub struct RefIndex {
    path: PathBuf,
    handle: Mutex<FaidxHandle>,
}

impl RefIndex {
    /// Load a reference FASTA, building `.fai` if absent.
    pub fn load(path: impl Into<PathBuf>) -> Result<Self, String> {
        ensure_default_htslib_log_level();
        let path = path.into();
        let handle = FaidxHandle::load(path.to_str().ok_or("non-UTF8 ref path")?)
            .ok_or_else(|| format!("fai_load failed for {}", path.display()))?;
        Ok(RefIndex {
            path,
            handle: Mutex::new(handle),
        })
    }

    pub fn path(&self) -> &std::path::Path {
        &self.path
    }

    pub fn nseq(&self) -> i32 {
        self.handle.lock().unwrap().nseq()
    }

    pub fn has_seq(&self, name: &str) -> bool {
        self.handle.lock().unwrap().has_seq(name)
    }

    pub fn seq_len(&self, name: &str) -> Option<i64> {
        self.handle.lock().unwrap().seq_len(name)
    }

    /// Fetch `[beg, end]` (**0-based, inclusive**) as owned bytes, plus strand.
    /// Returns the plus-strand sequence. `beg` is clamped to >= 0.
    pub fn fetch_0based(&self, chr: &str, beg: i64, end: i64) -> Result<Vec<u8>, String> {
        let beg = beg.max(0);
        self.handle.lock().unwrap().fetch(chr, beg, end)
    }

    /// Fetch `[start, end]` (**1-based, inclusive**) as owned bytes.
    /// Equivalent to `samtools faidx chr:start-end`.
    pub fn fetch_1based(&self, chr: &str, start: i64, end: i64) -> Result<Vec<u8>, String> {
        if start < 1 {
            return Err(format!("1-based start {} < 1", start));
        }
        self.fetch_0based(chr, start - 1, end - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ref() -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("consensus_rs_refidx_test");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // chr1: 50 A | 50 C | 50 G | 50 T  (200 bp)
        let seq = format!(
            "{}{}{}{}",
            "A".repeat(50),
            "C".repeat(50),
            "G".repeat(50),
            "T".repeat(50)
        );
        let fa = dir.join("ref.fa");
        std::fs::write(&fa, format!(">chr1\n{}\n", seq)).unwrap();
        fa
    }

    #[test]
    fn fetch_1based_matches_samtools_semantics() {
        let fa = make_ref();
        let r = RefIndex::load(&fa).unwrap();
        assert_eq!(r.nseq(), 1);
        assert_eq!(r.seq_len("chr1"), Some(200));

        // 1-based [1,100] == 0-based [0,99] -> 100 bp, all A then C
        let s = r.fetch_1based("chr1", 1, 100).unwrap();
        assert_eq!(s.len(), 100);
        assert!(s[..50].iter().all(|&b| b == b'A'));
        assert!(s[50..].iter().all(|&b| b == b'C'));

        // 1-based [101,200] -> 50 G + 50 T
        let s = r.fetch_1based("chr1", 101, 200).unwrap();
        assert_eq!(s.len(), 100);
        assert!(s[..50].iter().all(|&b| b == b'G'));
        assert!(s[50..].iter().all(|&b| b == b'T'));

        let _ = std::fs::remove_dir_all(fa.parent().unwrap());
    }
}
