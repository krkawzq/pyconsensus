//! consensus-rs — Rust rewrite of `bcftools consensus`.
//!
//! See `docs/design.md` and `docs/implementation_plan.md`.
//!
//! Layout (filled in across milestones M0–M5):
//!   - `htslib_ffi` : `extern "C"` bindings to htslib (faidx / vcf / bcf / regidx)
//!   - `ref_index`  : faidx wrapper, per-region fetch               (M1)
//!   - `vcf_store`  : one-shot VCF preprocessing + region query     (M1)
//!   - `iupac`      : iupac tables ported from bcftools.h           (M3)
//!   - `haplotype`  : `-H` / sample-mode allele selection           (M3)
//!   - `mask`/`chain`: `-m` / `-c`                                  (M4)
//!   - `apply`      : consensus apply state machine (Rust rewrite)  (M2–M4)
//!   - `engine`     : thread pool + channel                         (M5)
//!   - `py`         : PyO3 bindings (feature `python`)              (M5)

pub mod apply;
pub mod chain;
pub mod compiled;
pub mod engine;
pub mod haplotype;
pub mod htslib_ffi;
pub mod iupac;
pub mod logging;
pub mod mask;
pub mod planner;
pub mod ref_index;
pub mod stats;
pub mod vcf_store;

#[cfg(feature = "python")]
pub mod py;

#[cfg(test)]
mod tests {
    use super::htslib_ffi::FaidxHandle;
    use std::fs;

    /// M0 acceptance: load a ref, fetch a subsequence, confirm the 0-based
    /// **closed** interval semantics of `faidx_fetch_seq64`.
    #[test]
    fn faidx_fetch_subsequence_closed_interval() {
        let dir = std::env::temp_dir().join("consensus_rs_m0_test");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        // 200 bp: 50 A | 50 C | 50 G | 50 T
        let seq = format!(
            "{}{}{}{}",
            "A".repeat(50),
            "C".repeat(50),
            "G".repeat(50),
            "T".repeat(50)
        );
        let fa = dir.join("ref.fa");
        fs::write(&fa, format!(">chr1\n{}\n", seq)).unwrap();

        let fai = FaidxHandle::load(fa.to_str().unwrap()).expect("fai_load failed");

        assert_eq!(fai.nseq(), 1);
        assert!(fai.has_seq("chr1"));
        assert_eq!(fai.seq_len("chr1"), Some(200));

        // 0-based closed [0, 99] -> 100 bp
        let s = fai.fetch("chr1", 0, 99).unwrap();
        assert_eq!(s.len(), 100, "[0,99] must be 100 bp (closed interval)");
        // first 50 are A, next 50 are C
        assert!(s[..50].iter().all(|&b| b == b'A'));
        assert!(s[50..].iter().all(|&b| b == b'C'));

        // [0, 0] -> single base
        assert_eq!(
            fai.fetch("chr1", 0, 0).unwrap().len(),
            1,
            "[0,0] must be 1 bp"
        );

        // boundary [49, 50] -> A then C
        assert_eq!(fai.fetch("chr1", 49, 50).unwrap(), b"AC");

        // tail [150, 199] -> 50 T
        let t = fai.fetch("chr1", 150, 199).unwrap();
        assert_eq!(t.len(), 50);
        assert!(t.iter().all(|&b| b == b'T'));

        let _ = fs::remove_dir_all(&dir);
    }
}
