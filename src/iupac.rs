//! iupac — IUPAC nucleotide code tables, ported verbatim from
//! `bcftools/bcftools.h:107-136` (`iupac2bitmask` / `bitmask2iupac`).
//!
//! These are pure lookup tables with no htslib dependency.

use smallvec::SmallVec;

pub type IupacAlleleBuf = SmallVec<[u8; 64]>;

const A: i32 = 1;
const C: i32 = 2;
const G: i32 = 4;
const T: i32 = 8;

/// Map an IUPAC character to its base bitmask (A|C|G|T bits).
/// Returns -1 for symbolic / invalid characters. Case-insensitive.
pub fn iupac2bitmask(mut iupac: u8) -> i32 {
    if iupac >= b'a' {
        iupac -= 32;
    }
    match iupac {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        b'M' => A | C,
        b'R' => A | G,
        b'W' => A | T,
        b'S' => C | G,
        b'Y' => C | T,
        b'K' => G | T,
        b'V' => A | C | G,
        b'H' => A | C | T,
        b'D' => A | G | T,
        b'B' => C | G | T,
        b'N' => A | C | G | T,
        _ => -1,
    }
}

/// Inverse of `iupac2bitmask`. Returns 0 for an invalid bitmask.
pub fn bitmask2iupac(bitmask: i32) -> u8 {
    const IUPAC: [u8; 16] = [
        b'.', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B',
        b'N',
    ];
    if bitmask <= 0 || bitmask > 15 {
        return 0;
    }
    IUPAC[bitmask as usize]
}

/// `iupac_set_allele` (consensus.c:552-581, also 730-758 for the no-sample
/// `-I` path): combine a set of allele strings by OR-ing per-position bitmasks,
/// writing the resulting IUPAC bytes into the longest allele, and returning its
/// index.
///
/// `alleles` is REF + ALTs (indices 0..). `selected` marks which alleles
/// (by index) to mix in. Returns `(ialt, modified_alt_bytes)` where `ialt` is
/// the index of the allele that received the IUPAC bytes (the longest selected
/// one with index > 0), or a fallback if no ALT qualifies.
///
/// Mirrors the no-sample `output_iupac` block at consensus.c:730-758.
pub fn iupac_set_allele(alleles: &[&[u8]], selected: &[bool]) -> (Option<usize>, IupacAlleleBuf) {
    iupac_set_allele_impl(alleles, |i| selected.get(i).copied().unwrap_or(false))
}

pub fn iupac_set_all_alleles(alleles: &[&[u8]]) -> (Option<usize>, IupacAlleleBuf) {
    iupac_set_allele_impl(alleles, |_| true)
}

pub fn iupac_set_allele_mask(
    alleles: &[&[u8]],
    selected_mask: u64,
) -> (Option<usize>, IupacAlleleBuf) {
    iupac_set_allele_impl(alleles, |i| i < 64 && (selected_mask & (1u64 << i)) != 0)
}

fn iupac_set_allele_impl(
    alleles: &[&[u8]],
    mut is_selected: impl FnMut(usize) -> bool,
) -> (Option<usize>, IupacAlleleBuf) {
    let mut selected = SmallVec::<[bool; 8]>::with_capacity(alleles.len());
    for i in 0..alleles.len() {
        selected.push(is_selected(i));
    }
    if let Some(out) = iupac_set_single_base_fast(alleles, &selected) {
        return out;
    }

    let mut max_len: usize = 0;
    let mut alt_len: usize = 0;
    let mut ialt: i32 = -1;
    let mut fallback_alt: i32 = -1;
    let mut bitmask: SmallVec<[i32; 64]> = SmallVec::new();

    for (i, al) in alleles.iter().enumerate() {
        if !selected[i] {
            continue;
        }
        if fallback_alt <= 0 {
            fallback_alt = i as i32;
        }
        // skip symbolic / invalid
        let mut j = 0;
        while j < al.len() {
            if iupac2bitmask(al[j]) < 0 {
                break;
            }
            j += 1;
        }
        if j < al.len() {
            continue;
        }
        let l = al.len();
        if l > max_len {
            bitmask.resize(l, 0);
            max_len = l;
        }
        if i > 0 && l > alt_len {
            alt_len = l;
            ialt = i as i32;
        }
        for j in 0..l {
            bitmask[j] |= iupac2bitmask(al[j]);
        }
    }

    if alt_len > 0 {
        // write IUPAC bytes into a copy of the longest allele
        let mut out = IupacAlleleBuf::from_slice(alleles[ialt as usize]);
        for j in 0..alt_len {
            out[j] = bitmask2iupac(bitmask[j]);
        }
        return (Some(ialt as usize), out);
    }
    if fallback_alt >= 0 {
        let i = fallback_alt as usize;
        return (Some(i), IupacAlleleBuf::from_slice(alleles[i]));
    }
    (None, IupacAlleleBuf::new())
}

fn iupac_set_single_base_fast(
    alleles: &[&[u8]],
    selected: &[bool],
) -> Option<(Option<usize>, IupacAlleleBuf)> {
    let mut bitmask = 0;
    let mut ialt: Option<usize> = None;
    let mut fallback_alt: Option<usize> = None;
    let mut any_selected = false;

    for (i, al) in alleles.iter().enumerate() {
        if !selected.get(i).copied().unwrap_or(false) {
            continue;
        }
        any_selected = true;
        if fallback_alt.is_none() {
            fallback_alt = Some(i);
        }
        if al.len() != 1 {
            return None;
        }
        let m = iupac2bitmask(al[0]);
        if m < 0 {
            return None;
        }
        bitmask |= m;
        if i > 0 && ialt.is_none() {
            ialt = Some(i);
        }
    }

    if !any_selected {
        return Some((None, IupacAlleleBuf::new()));
    }
    if let Some(i) = ialt {
        let iupac = bitmask2iupac(bitmask);
        let mut out = IupacAlleleBuf::new();
        if iupac == 0 {
            return None;
        }
        out.push(iupac);
        return Some((Some(i), out));
    }
    fallback_alt.map(|i| (Some(i), IupacAlleleBuf::from_slice(alleles[i])))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_single_bases() {
        for &b in b"ACGT".iter() {
            assert_eq!(bitmask2iupac(iupac2bitmask(b)), b);
        }
    }

    #[test]
    fn mixed_iupac() {
        // A|C -> M
        assert_eq!(
            bitmask2iupac(iupac2bitmask(b'A') | iupac2bitmask(b'C')),
            b'M'
        );
        // A|C|G|T -> N
        assert_eq!(
            bitmask2iupac(
                iupac2bitmask(b'A')
                    | iupac2bitmask(b'C')
                    | iupac2bitmask(b'G')
                    | iupac2bitmask(b'T')
            ),
            b'N'
        );
        // lowercase
        assert_eq!(iupac2bitmask(b'g'), iupac2bitmask(b'G'));
    }

    #[test]
    fn symbolic_returns_neg1() {
        assert_eq!(iupac2bitmask(b'<'), -1);
        assert_eq!(iupac2bitmask(b'*'), -1);
    }

    #[test]
    fn set_allele_mixed_ref_alt() {
        // REF=A, ALT=C, both selected -> IUPAC M
        let alleles: Vec<&[u8]> = vec![b"A", b"C"];
        let sel = vec![true, true];
        let (ialt, out) = iupac_set_allele(&alleles, &sel);
        assert_eq!(ialt, Some(1));
        assert_eq!(out.as_slice(), b"M");
    }

    #[test]
    fn set_allele_skips_symbolic() {
        // REF=A, ALT=<DEL> (symbolic, skipped), ALT2=C -> mix A and C
        let alleles: Vec<&[u8]> = vec![b"A", b"<DEL>", b"C"];
        let sel = vec![true, true, true];
        let (ialt, out) = iupac_set_allele(&alleles, &sel);
        assert_eq!(ialt, Some(2));
        assert_eq!(out.as_slice(), b"M");
    }
}
