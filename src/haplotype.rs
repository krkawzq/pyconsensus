//! haplotype — `-H` parsing, sample-mode classification, and allele selection.
//!
//! Ports `main_consensus` `-H` parsing (consensus.c:1310-1328), the sample /
//! haplotype mode decision in `init_data` (consensus.c:250-269), and the
//! GT-driven allele selection in `apply_variant` (consensus.c:602-758).

use crate::iupac::{
    iupac_set_all_alleles, iupac_set_allele, iupac_set_allele_mask, IupacAlleleBuf,
};
use crate::vcf_store::{GtAllele, VcfRecord};
use smallvec::SmallVec;

// PICK_* flags (consensus.c:51-55)
pub const PICK_REF: u8 = 1;
pub const PICK_ALT: u8 = 2;
pub const PICK_LONG: u8 = 4;
pub const PICK_SHORT: u8 = 8;
pub const PICK_IUPAC: u8 = 16;

/// Parsed `-H` spec: a PICK_* bitmask plus an optional haplotype index.
#[derive(Clone, Debug, Default)]
pub struct HaplotypeSpec {
    pub pick: u8,
    pub haplotype: Option<u32>,
}

impl HaplotypeSpec {
    /// Parse a cli `-H` string. Returns None on unrecognised input.
    /// Mirrors consensus.c:1310-1328.
    pub fn parse(s: &str) -> Option<Self> {
        let upper = s.to_ascii_uppercase();
        let mut spec = HaplotypeSpec::default();
        match upper.as_str() {
            "R" => spec.pick |= PICK_REF,
            "A" => spec.pick |= PICK_ALT,
            "L" => spec.pick |= PICK_LONG | PICK_REF,
            "S" => spec.pick |= PICK_SHORT | PICK_REF,
            "LR" => spec.pick |= PICK_LONG | PICK_REF,
            "LA" => spec.pick |= PICK_LONG | PICK_ALT,
            "SR" => spec.pick |= PICK_SHORT | PICK_REF,
            "SA" => spec.pick |= PICK_SHORT | PICK_ALT,
            "I" => spec.pick |= PICK_IUPAC,
            "1PIU" => {
                spec.pick |= PICK_IUPAC;
                spec.haplotype = Some(1);
            }
            "2PIU" => {
                spec.pick |= PICK_IUPAC;
                spec.haplotype = Some(2);
            }
            _ => {
                // "<N>" or "<N>pIu"
                let (num_part, rest) = match upper.find(|c: char| !c.is_ascii_digit()) {
                    Some(i) => (&upper[..i], &upper[i..]),
                    None => (upper.as_str(), ""),
                };
                if num_part.is_empty() {
                    return None;
                }
                let n: u32 = num_part.parse().ok()?;
                if n == 0 {
                    return None;
                }
                spec.haplotype = Some(n);
                if rest.is_empty() {
                    // plain number -> use_hap, no IUPAC
                } else if rest == "PIU" {
                    spec.pick |= PICK_IUPAC;
                } else {
                    return None;
                }
            }
        }
        Some(spec)
    }
}

/// How alleles are chosen for a record, mirroring init_data's classification.
#[derive(Clone, Debug)]
pub enum SampleMode {
    ApplyAllAlt,
    IupacAllSamples { samples: Vec<i32> },
    SingleSample { idx: i32, spec: HaplotypeSpec },
    IupacFromRefAlt,
}

impl Default for SampleMode {
    fn default() -> Self {
        SampleMode::ApplyAllAlt
    }
}

/// Result of selecting an allele for a record.
/// `ialt = None` means "skip this record" (bcftools `return`).
/// `ialt = Some(0)` means REF (freeze). `ialt = Some(i)` is the ALT index.
/// When IUPAC mixing produces modified bytes, they're carried in `alt_override`
/// and the apply step uses them instead of `rec.alleles[ialt]`.
#[derive(Clone, Debug)]
pub struct AlleleSelection {
    pub ialt: Option<i32>,
    pub alt_override: Option<IupacAlleleBuf>,
}

/// Select the allele for `rec` under `mode`. `missing_allele` controls the
/// missing-GT path (Some → emit char, None → skip). Mirrors apply_variant
/// 602-758.
#[inline]
pub fn select_allele(
    rec: &VcfRecord,
    mode: &SampleMode,
    missing_allele: Option<u8>,
) -> AlleleSelection {
    match mode {
        SampleMode::ApplyAllAlt => AlleleSelection {
            ialt: Some(1),
            alt_override: None,
        },
        SampleMode::IupacFromRefAlt if rec.alleles.len() > 1 => {
            let alleles = allele_slices(rec);
            let (ialt, out) = iupac_set_all_alleles(&alleles);
            AlleleSelection {
                ialt: ialt.map(|i| i as i32),
                alt_override: if out.is_empty() { None } else { Some(out) },
            }
        }
        SampleMode::IupacFromRefAlt => AlleleSelection {
            ialt: None,
            alt_override: None,
        },
        SampleMode::IupacAllSamples { samples } => {
            // consensus.c:603-619: accumulate GT alleles across samples, IUPAC-mix.
            if rec.gt.is_empty() {
                return AlleleSelection {
                    ialt: None,
                    alt_override: None,
                };
            }
            let n_allele = rec.alleles.len();
            if n_allele <= 64 {
                let mut selected_mask = 0u64;
                for &s in samples {
                    let s = s as usize;
                    if s >= rec.gt.len() {
                        continue;
                    }
                    for a in &rec.gt[s] {
                        if let Some(idx) = a.allele {
                            let idx = idx as usize;
                            if idx < n_allele {
                                selected_mask |= 1u64 << idx;
                            }
                        }
                    }
                }
                if selected_mask == 0 {
                    if missing_allele.is_none() {
                        return AlleleSelection {
                            ialt: None,
                            alt_override: None,
                        };
                    }
                    return AlleleSelection {
                        ialt: Some(-1),
                        alt_override: None,
                    };
                }
                let alleles = allele_slices(rec);
                let (ialt, out) = iupac_set_allele_mask(&alleles, selected_mask);
                return AlleleSelection {
                    ialt: ialt.map(|i| i as i32),
                    alt_override: if out.is_empty() { None } else { Some(out) },
                };
            }
            let mut selected: SmallVec<[bool; 8]> = SmallVec::new();
            selected.resize(n_allele, false);
            let mut is_set = false;
            for &s in samples {
                let s = s as usize;
                if s >= rec.gt.len() {
                    continue;
                }
                for a in &rec.gt[s] {
                    match a.allele {
                        None => continue, // missing
                        Some(idx) if (idx as usize) < n_allele => {
                            selected[idx as usize] = true;
                            is_set = true;
                        }
                        _ => {}
                    }
                }
            }
            if !is_set {
                if missing_allele.is_none() {
                    return AlleleSelection {
                        ialt: None,
                        alt_override: None,
                    };
                }
                return AlleleSelection {
                    ialt: Some(-1),
                    alt_override: None,
                };
            }
            let alleles = allele_slices(rec);
            let (ialt, out) = iupac_set_allele(&alleles, &selected);
            AlleleSelection {
                ialt: ialt.map(|i| i as i32),
                alt_override: if out.is_empty() { None } else { Some(out) },
            }
        }
        SampleMode::SingleSample { idx, spec } => {
            select_single_sample(rec, *idx, spec, missing_allele)
        }
    }
}

/// Single-sample selection (consensus.c:622-758).
fn select_single_sample(
    rec: &VcfRecord,
    idx: i32,
    spec: &HaplotypeSpec,
    missing_allele: Option<u8>,
) -> AlleleSelection {
    if let Some(selection) = select_biallelic_hap_fast(rec, idx, spec, missing_allele) {
        return selection;
    }

    if rec.gt.is_empty() {
        return AlleleSelection {
            ialt: None,
            alt_override: None,
        };
    }
    let gt = match rec.gt.get(idx as usize) {
        Some(g) => g,
        None => {
            return AlleleSelection {
                ialt: None,
                alt_override: None,
            }
        }
    };
    let ploidy = gt.len();

    // Decide action (consensus.c:631-640)
    enum Action {
        UseHap,
        UseIupac,
        PickOne,
    }
    let action = if spec.pick & PICK_IUPAC != 0 {
        if spec.haplotype.is_none() {
            Action::UseIupac
        } else if !is_phased(gt, 0) && !is_phased_last(gt) {
            Action::UseIupac
        } else {
            Action::UseHap
        }
    } else {
        // (output_iupac without -H is IupacFromRefAlt, not this branch)
        if spec.haplotype.is_none() {
            Action::PickOne
        } else {
            Action::UseHap
        }
    };

    match action {
        Action::UseHap => {
            let hap = spec.haplotype.unwrap_or(1) as usize;
            if hap > ploidy {
                // missing on either end?
                let missing_end = gt.last().map(|a| a.allele.is_none()).unwrap_or(true)
                    || gt.first().map(|a| a.allele.is_none()).unwrap_or(true);
                if missing_end {
                    if missing_allele.is_none() {
                        return AlleleSelection {
                            ialt: None,
                            alt_override: None,
                        };
                    }
                    return AlleleSelection {
                        ialt: Some(-1),
                        alt_override: None,
                    };
                }
                // warn + skip (bcftools prints once)
                return AlleleSelection {
                    ialt: None,
                    alt_override: None,
                };
            }
            let a = &gt[hap - 1];
            if a.allele.is_none() {
                if missing_allele.is_none() {
                    return AlleleSelection {
                        ialt: None,
                        alt_override: None,
                    };
                }
                return AlleleSelection {
                    ialt: Some(-1),
                    alt_override: None,
                };
            }
            AlleleSelection {
                ialt: a.allele.map(|i| i as i32),
                alt_override: None,
            }
        }
        Action::UseIupac => {
            // IUPAC-mix this single sample's GT (consensus.c:671-678)
            let n_allele = rec.alleles.len();
            if n_allele <= 64 {
                let mut selected_mask = 0u64;
                for a in gt {
                    if let Some(idx) = a.allele {
                        let idx = idx as usize;
                        if idx < n_allele {
                            selected_mask |= 1u64 << idx;
                        }
                    }
                }
                if selected_mask == 0 {
                    if missing_allele.is_none() {
                        return AlleleSelection {
                            ialt: None,
                            alt_override: None,
                        };
                    }
                    return AlleleSelection {
                        ialt: Some(-1),
                        alt_override: None,
                    };
                }
                let alleles = allele_slices(rec);
                let (ialt, out) = iupac_set_allele_mask(&alleles, selected_mask);
                return AlleleSelection {
                    ialt: ialt.map(|i| i as i32),
                    alt_override: if out.is_empty() { None } else { Some(out) },
                };
            }
            let mut selected: SmallVec<[bool; 8]> = SmallVec::new();
            selected.resize(n_allele, false);
            let mut is_set = false;
            for a in gt {
                if a.allele.is_none() {
                    continue;
                }
                if let Some(idx) = a.allele {
                    if (idx as usize) < n_allele {
                        selected[idx as usize] = true;
                        is_set = true;
                    }
                }
            }
            if !is_set {
                if missing_allele.is_none() {
                    return AlleleSelection {
                        ialt: None,
                        alt_override: None,
                    };
                }
                return AlleleSelection {
                    ialt: Some(-1),
                    alt_override: None,
                };
            }
            let alleles = allele_slices(rec);
            let (ialt, out) = iupac_set_allele(&alleles, &selected);
            AlleleSelection {
                ialt: ialt.map(|i| i as i32),
                alt_override: if out.is_empty() { None } else { Some(out) },
            }
        }
        Action::PickOne => {
            // consensus.c:679-721: hom → that allele; het → PICK_LONG/SHORT/REF/ALT
            let mut ialt: i32 = 0;
            let mut is_hom = true;
            let mut first: Option<i32> = None;
            for a in gt {
                let cur = match a.allele {
                    Some(v) => v as i32,
                    None => {
                        if missing_allele.is_none() {
                            return AlleleSelection {
                                ialt: None,
                                alt_override: None,
                            };
                        }
                        return AlleleSelection {
                            ialt: Some(-1),
                            alt_override: None,
                        };
                    }
                };
                if let Some(first) = first {
                    if cur != first {
                        is_hom = false;
                        break;
                    }
                } else {
                    first = Some(cur);
                }
                ialt = cur;
            }
            if is_hom {
                return AlleleSelection {
                    ialt: Some(ialt),
                    alt_override: None,
                };
            }
            // het: apply PICK_LONG / PICK_SHORT / PICK_REF / PICK_ALT
            let mut prev_len: i64 = 0;
            let mut chosen: i32 = 0;
            for a in gt {
                let jalt = match a.allele {
                    Some(x) => x as i32,
                    None => continue,
                };
                if rec.alleles.len() <= jalt as usize {
                    continue;
                }
                let len = if jalt == 0 {
                    rec.rlen as i64
                } else {
                    rec.alleles[jalt as usize].len() as i64
                };
                if spec.pick & (PICK_LONG | PICK_SHORT) != 0 {
                    if prev_len == 0 {
                        chosen = jalt;
                        prev_len = len;
                    } else if len == prev_len {
                        if spec.pick & PICK_REF != 0 && jalt == 0 {
                            chosen = jalt;
                            prev_len = len;
                        } else if spec.pick & PICK_ALT != 0 && chosen == 0 {
                            chosen = jalt;
                            prev_len = len;
                        }
                    } else if spec.pick & PICK_LONG != 0 && len > prev_len {
                        chosen = jalt;
                        prev_len = len;
                    } else if spec.pick & PICK_SHORT != 0 && len < prev_len {
                        chosen = jalt;
                        prev_len = len;
                    }
                } else {
                    if spec.pick & PICK_REF != 0 && jalt == 0 {
                        chosen = jalt;
                    } else if spec.pick & PICK_ALT != 0 && chosen == 0 {
                        chosen = jalt;
                    }
                }
            }
            AlleleSelection {
                ialt: Some(chosen),
                alt_override: None,
            }
        }
    }
}

fn select_biallelic_hap_fast(
    rec: &VcfRecord,
    idx: i32,
    spec: &HaplotypeSpec,
    missing_allele: Option<u8>,
) -> Option<AlleleSelection> {
    if idx < 0 {
        return None;
    }
    let hap = spec.haplotype? as usize;
    if hap == 0 || hap > 2 {
        return None;
    }
    let bits = rec.gt_bits.as_ref()?;
    match bits.allele_for_hap(idx as usize, hap)? {
        Some(ialt) => Some(AlleleSelection {
            ialt: Some(ialt),
            alt_override: None,
        }),
        None if missing_allele.is_some() => Some(AlleleSelection {
            ialt: Some(-1),
            alt_override: None,
        }),
        None => Some(AlleleSelection {
            ialt: None,
            alt_override: None,
        }),
    }
}

fn allele_slices(rec: &VcfRecord) -> SmallVec<[&[u8]; 8]> {
    let mut alleles = SmallVec::with_capacity(rec.alleles.len());
    for allele in &rec.alleles {
        alleles.push(&allele[..]);
    }
    alleles
}

fn is_phased(gt: &[GtAllele], i: usize) -> bool {
    gt.get(i).map(|a| a.phased).unwrap_or(false)
}

fn is_phased_last(gt: &[GtAllele]) -> bool {
    gt.last().map(|a| a.phased).unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hap_specs() {
        assert_eq!(HaplotypeSpec::parse("R").unwrap().pick & PICK_REF, PICK_REF);
        assert_eq!(HaplotypeSpec::parse("A").unwrap().pick & PICK_ALT, PICK_ALT);
        let lr = HaplotypeSpec::parse("LR").unwrap();
        assert!(lr.pick & PICK_LONG != 0 && lr.pick & PICK_REF != 0);
        let i1 = HaplotypeSpec::parse("1pIu").unwrap();
        assert!(i1.pick & PICK_IUPAC != 0);
        assert_eq!(i1.haplotype, Some(1));
        let n3 = HaplotypeSpec::parse("3").unwrap();
        assert_eq!(n3.haplotype, Some(3));
        assert_eq!(n3.pick, 0);
        let n3piu = HaplotypeSpec::parse("3pIu").unwrap();
        assert!(n3piu.pick & PICK_IUPAC != 0);
        assert_eq!(n3piu.haplotype, Some(3));
        assert!(HaplotypeSpec::parse("garbage").is_none());
        assert!(HaplotypeSpec::parse("0").is_none());
    }

    #[test]
    fn apply_all_alt_mode() {
        use smallvec::SmallVec;
        // irrelevant record fields for ApplyAllAlt
        let rec = VcfRecord {
            pos: 0,
            rlen: 1,
            rid: 0,
            alleles: vec![SmallVec::from_slice(b"A"), SmallVec::from_slice(b"G")],
            gt: vec![],
            gt_bits: None,
            var_type: 1,
            compiled: crate::compiled::CompiledRecord::from_alleles(
                1,
                &[SmallVec::from_slice(b"A"), SmallVec::from_slice(b"G")],
            ),
        };
        let sel = select_allele(&rec, &SampleMode::ApplyAllAlt, None);
        assert_eq!(sel.ialt, Some(1));
        assert!(sel.alt_override.is_none());
    }

    #[test]
    fn single_sample_haplotype_uses_biallelic_bitset() {
        use smallvec::SmallVec;

        let mut s0: SmallVec<[GtAllele; 2]> = SmallVec::new();
        s0.push(GtAllele {
            allele: Some(0),
            phased: true,
            raw: 0,
        });
        s0.push(GtAllele {
            allele: Some(1),
            phased: true,
            raw: 0,
        });
        let mut s1: SmallVec<[GtAllele; 2]> = SmallVec::new();
        s1.push(GtAllele {
            allele: Some(1),
            phased: true,
            raw: 0,
        });
        s1.push(GtAllele {
            allele: Some(0),
            phased: true,
            raw: 0,
        });
        let gt = vec![s0, s1];
        let gt_bits = crate::vcf_store::BiallelicPhasedGtBits::from_gt(2, 2, &gt);
        let rec = VcfRecord {
            pos: 0,
            rlen: 1,
            rid: 0,
            alleles: vec![SmallVec::from_slice(b"A"), SmallVec::from_slice(b"G")],
            gt,
            gt_bits,
            var_type: 1,
            compiled: crate::compiled::CompiledRecord::from_alleles(
                1,
                &[SmallVec::from_slice(b"A"), SmallVec::from_slice(b"G")],
            ),
        };

        let hap1 = SampleMode::SingleSample {
            idx: 0,
            spec: HaplotypeSpec::parse("1").unwrap(),
        };
        let hap2 = SampleMode::SingleSample {
            idx: 0,
            spec: HaplotypeSpec::parse("2").unwrap(),
        };

        assert_eq!(select_allele(&rec, &hap1, None).ialt, Some(0));
        assert_eq!(select_allele(&rec, &hap2, None).ialt, Some(1));
    }
}
