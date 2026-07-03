//! mask — compiled interval store for `-m FILE` + `--mask-with CHAR|uc|lc`.
//!
//! The hot path keeps BED intervals in Rust-owned vectors keyed by contig. This
//! avoids per-region htslib regidx iterators, CString creation, and mutable FFI
//! iterator state.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Interval {
    beg: i64,
    end: i64,
}

/// One compiled mask file + replacement mode.
pub struct Mask {
    by_chr: HashMap<String, Vec<Interval>>,
    pub with: MaskWith,
    path: PathBuf,
}

impl Mask {
    /// Load a BED-like mask file. Coordinates are BED 0-based half-open
    /// `[start,end)` and are stored internally as inclusive intervals.
    pub fn load(path: impl Into<PathBuf>, with: MaskWith) -> Result<Self, String> {
        let path = path.into();
        let f = File::open(&path).map_err(|e| format!("open mask {}: {}", path.display(), e))?;
        let mut by_chr: HashMap<String, Vec<Interval>> = HashMap::new();
        for (line_no, line) in BufReader::new(f).lines().enumerate() {
            let line = line.map_err(|e| format!("read mask {}: {}", path.display(), e))?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut fields = line.split_whitespace();
            let chr = match fields.next() {
                Some(chr) => chr,
                None => continue,
            };
            let beg: i64 = fields
                .next()
                .ok_or_else(|| format!("mask {}:{} missing start", path.display(), line_no + 1))?
                .parse()
                .map_err(|_| format!("mask {}:{} invalid start", path.display(), line_no + 1))?;
            let end_excl: i64 = fields
                .next()
                .ok_or_else(|| format!("mask {}:{} missing end", path.display(), line_no + 1))?
                .parse()
                .map_err(|_| format!("mask {}:{} invalid end", path.display(), line_no + 1))?;
            if beg < 0 || end_excl < beg {
                return Err(format!(
                    "mask {}:{} invalid interval",
                    path.display(),
                    line_no + 1
                ));
            }
            if end_excl == beg {
                continue;
            }
            by_chr.entry(chr.to_string()).or_default().push(Interval {
                beg,
                end: end_excl - 1,
            });
        }
        for intervals in by_chr.values_mut() {
            intervals.sort_by_key(|iv| (iv.beg, iv.end));
            let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
            for iv in intervals.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if iv.beg <= last.end + 1 {
                        if iv.end > last.end {
                            last.end = iv.end;
                        }
                        continue;
                    }
                }
                merged.push(iv);
            }
            *intervals = merged;
        }
        Ok(Mask { by_chr, with, path })
    }

    pub fn path(&self) -> &std::path::Path {
        &self.path
    }

    /// Does any mask region overlap `[start, end]` (0-based, inclusive)?
    pub fn overlaps(&self, chr: &str, start: i64, end: i64) -> bool {
        if end < start {
            return false;
        }
        let Some(intervals) = self.by_chr.get(chr) else {
            return false;
        };
        let first_after_end = intervals.partition_point(|iv| iv.beg <= end);
        if first_after_end == 0 {
            return false;
        }
        intervals[first_after_end - 1].end >= start
    }

    /// Apply the mask to `buf`, which spans `[ori_pos, ori_pos+len-1]`
    /// (0-based, inclusive).
    pub fn apply_to_buf(&self, chr: &str, buf: &mut [u8], ori_pos: i64) {
        if buf.is_empty() {
            return;
        }
        let Some(intervals) = self.by_chr.get(chr) else {
            return;
        };
        let len = buf.len() as i64;
        let start = ori_pos;
        let end = ori_pos + len - 1;
        let mut i = intervals.partition_point(|iv| iv.end < start);
        while i < intervals.len() {
            let iv = intervals[i];
            if iv.beg > end {
                break;
            }
            let s = (iv.beg.max(start) - start) as usize;
            let e = (iv.end.min(end) - start) as usize;
            match self.with {
                MaskWith::Char(ch) => buf[s..=e].fill(ch),
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
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_mask(name: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("consensus_rs_mask_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("mask.bed");
        std::fs::write(&path, body).unwrap();
        path
    }

    #[test]
    fn compiled_mask_overlaps_and_applies_bed_intervals() {
        let path = write_mask("basic", "# comment\nchr1\t2\t5\nchr1\t5\t7\nchr1\t10\t12\n");
        let mask = Mask::load(&path, MaskWith::Char(b'N')).unwrap();

        assert!(mask.overlaps("chr1", 1, 2));
        assert!(mask.overlaps("chr1", 6, 6));
        assert!(!mask.overlaps("chr1", 7, 9));
        assert!(!mask.overlaps("chr2", 2, 4));

        let mut buf = b"abcdefghijklmn".to_vec();
        mask.apply_to_buf("chr1", &mut buf, 0);
        assert_eq!(buf, b"abNNNNNhijNNmn");
    }

    #[test]
    fn compiled_mask_case_modes() {
        let path = write_mask("case", "chr1\t1\t4\n");
        let mut buf = b"aCgTa".to_vec();
        Mask::load(&path, MaskWith::Uc)
            .unwrap()
            .apply_to_buf("chr1", &mut buf, 0);
        assert_eq!(buf, b"aCGTa");

        let mut buf = b"aCgTa".to_vec();
        Mask::load(&path, MaskWith::Lc)
            .unwrap()
            .apply_to_buf("chr1", &mut buf, 0);
        assert_eq!(buf, b"acgta");
    }
}
