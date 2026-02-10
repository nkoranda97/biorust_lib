use crate::alphabets::dna;
use crate::error::{BioError, BioResult};
use crate::seq::bytes::{self, IntoNeedle, Needle};
use crate::seq::protein::ProteinSeq;
use crate::seq::rna::RnaSeq;
use crate::seq::traits::SeqBytes;
use crate::seq::{best_frame_index, TranslationFrame};

use std::sync::LazyLock;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DnaSeq {
    bytes: Vec<u8>,
}

pub trait ReverseComplement: Sized {
    fn reverse_complement(&self) -> Self;
}

impl DnaSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        let alphabet = dna::iupac_alphabet();
        for (pos, &b) in bytes.iter().enumerate() {
            if !alphabet.symbols.contains(b as usize) {
                return Err(BioError::InvalidChar { ch: b as char, pos });
            }
        }
        Ok(Self { bytes })
    }

    #[inline]
    pub(crate) fn from_bytes_unchecked(bytes: Vec<u8>) -> Self {
        Self { bytes }
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub fn gc_content(&self) -> f64 {
        if self.bytes.is_empty() {
            return 0.0;
        }
        let gc = self
            .bytes
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
            .count();
        gc as f64 / self.bytes.len() as f64
    }

    pub fn reverse(&self) -> Self {
        let mut out = self.bytes.clone();
        out.reverse();
        Self { bytes: out }
    }

    pub fn reverse_complement(&self) -> Self {
        let out = dna::reverse_complement(self.as_bytes());
        Self { bytes: out }
    }

    pub fn complement(&self) -> Self {
        let mut out = Vec::with_capacity(self.bytes.len());
        for &base in self.as_bytes() {
            out.push(dna::complement(base));
        }
        Self { bytes: out }
    }

    pub fn transcribe(&self) -> RnaSeq {
        let mut out = self.bytes.clone();
        for b in &mut out {
            if *b == b'T' {
                *b = b'U';
            } else if *b == b't' {
                *b = b'u';
            }
        }
        RnaSeq::from_bytes_unchecked(out)
    }

    pub fn translate(&self) -> BioResult<ProteinSeq> {
        let bytes = self.as_bytes();
        if bytes.len() % 3 != 0 {
            return Err(BioError::TranslationError {
                msg: format!(
                    "sequence length {} is not a multiple of 3 ({} trailing bases would be lost)",
                    bytes.len(),
                    bytes.len() % 3
                ),
            });
        }
        Ok(translate_bytes(bytes, &BASE_INDEX))
    }

    pub fn translate_frame(&self, frame: TranslationFrame) -> BioResult<ProteinSeq> {
        match frame {
            TranslationFrame::One => {
                let bytes = self.as_bytes();
                let len = bytes.len() / 3 * 3;
                Ok(translate_bytes(&bytes[..len], &BASE_INDEX))
            }
            TranslationFrame::Two => {
                let bytes = self.as_bytes();
                if bytes.len() < 2 {
                    return Ok(ProteinSeq::from_bytes_unchecked(Vec::new()));
                }
                let slice = &bytes[1..];
                let len = slice.len() / 3 * 3;
                Ok(translate_bytes(&slice[..len], &BASE_INDEX))
            }
            TranslationFrame::Three => {
                let bytes = self.as_bytes();
                if bytes.len() < 3 {
                    return Ok(ProteinSeq::from_bytes_unchecked(Vec::new()));
                }
                let slice = &bytes[2..];
                let len = slice.len() / 3 * 3;
                Ok(translate_bytes(&slice[..len], &BASE_INDEX))
            }
            TranslationFrame::Auto => {
                let bytes = self.as_bytes();
                let mut candidates: [Vec<u8>; 3] = [Vec::new(), Vec::new(), Vec::new()];
                for offset in 0..3 {
                    if bytes.len() > offset {
                        let slice = &bytes[offset..];
                        let len = slice.len() / 3 * 3;
                        candidates[offset] = translate_to_vec(&slice[..len], &BASE_INDEX);
                    }
                }
                let idx = best_frame_index([
                    &candidates[0],
                    &candidates[1],
                    &candidates[2],
                ]);
                Ok(ProteinSeq::from_bytes_unchecked(
                    std::mem::take(&mut candidates[idx]),
                ))
            }
        }
    }

    pub fn count<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoNeedle<'a>,
    {
        let hay = self.as_bytes();
        let needle = sub.into_needle()?;

        Ok(bytes::count(hay, needle))
    }

    pub fn count_overlap<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoNeedle<'a>,
    {
        let hay = self.as_bytes();
        let needle = sub.into_needle()?;

        Ok(bytes::count_overlap(hay, needle))
    }
    pub fn contains<'a, N>(&'a self, sub: N) -> BioResult<bool>
    where
        N: IntoNeedle<'a>,
    {
        let hay = self.as_bytes();
        let needle = sub.into_needle()?;

        Ok(bytes::contains(hay, needle))
    }

    pub fn find<'a, N>(&'a self, sub: N, start: usize, end: usize) -> BioResult<Option<usize>>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::find(self.as_bytes(), needle, start, end))
    }

    pub fn rfind<'a, N>(&'a self, sub: N, start: usize, end: usize) -> BioResult<Option<usize>>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::rfind(self.as_bytes(), needle, start, end))
    }
}

impl SeqBytes for DnaSeq {
    fn as_bytes(&self) -> &[u8] {
        DnaSeq::as_bytes(self)
    }

    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self> {
        DnaSeq::new(bytes)
    }
}

impl ReverseComplement for DnaSeq {
    fn reverse_complement(&self) -> Self {
        DnaSeq::reverse_complement(self)
    }
}

impl<'a> IntoNeedle<'a> for &'a DnaSeq {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_bytes()))
    }
}

fn translate_to_vec(bytes: &[u8], base_index: &[u8; 256]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bytes.len() / 3);
    for codon in bytes.chunks_exact(3) {
        let i1 = base_index[codon[0] as usize];
        let i2 = base_index[codon[1] as usize];
        let i3 = base_index[codon[2] as usize];
        let aa = if i1 < 4 && i2 < 4 && i3 < 4 {
            let idx = ((i1 as usize) << 4) | ((i2 as usize) << 2) | (i3 as usize);
            CODON_TABLE[idx]
        } else {
            b'X'
        };
        out.push(aa);
    }
    out
}

fn translate_bytes(bytes: &[u8], base_index: &[u8; 256]) -> ProteinSeq {
    ProteinSeq::from_bytes_unchecked(translate_to_vec(bytes, base_index))
}

static BASE_INDEX: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut map = [255u8; 256];
    map[b'A' as usize] = 0;
    map[b'C' as usize] = 1;
    map[b'G' as usize] = 2;
    map[b'T' as usize] = 3;
    map[b'a' as usize] = 0;
    map[b'c' as usize] = 1;
    map[b'g' as usize] = 2;
    map[b't' as usize] = 3;
    map
});

const CODON_TABLE: [u8; 64] = *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_byte() {
        let s = DnaSeq::new(b"ACGTACGT".to_vec()).unwrap();
        assert_eq!(s.count(b'A').unwrap(), 2);
        assert_eq!(s.count(65u8).unwrap(), 2); // 'A'
        assert_eq!(s.count(65usize).unwrap(), 2);
        assert!(s.count(999usize).is_err());
    }

    #[test]
    fn count_subslice_nonoverlapping() {
        let s = DnaSeq::new(b"AAAAA".to_vec()).unwrap();
        assert_eq!(s.count(b"AA").unwrap(), 2); // non-overlapping: AA|AA|A
    }

    #[test]
    fn count_empty_subslice() {
        let s = DnaSeq::new(b"ACGT".to_vec()).unwrap();
        assert_eq!(s.count(b"").unwrap(), 5);
        assert_eq!(s.count("").unwrap(), 5);
    }

    #[test]
    fn count_from_other_dnaseq() {
        let s = DnaSeq::new(b"ACGTACGT".to_vec()).unwrap();
        let sub = DnaSeq::new(b"AC".to_vec()).unwrap();
        assert_eq!(s.count(&sub).unwrap(), 2);
    }

    #[test]
    fn contains_tests() {
        let s = DnaSeq::new(b"ACGTACGT".to_vec()).unwrap();

        assert!(s.contains(b"A").unwrap());
        assert!(s.contains(b"CG").unwrap());
        assert!(!s.contains(b"TTT").unwrap());
        assert!(s.contains(b"").unwrap());
    }

    #[test]
    fn find_basic() {
        let s = DnaSeq::new(b"ACGTACGT".to_vec()).unwrap();

        assert_eq!(s.find(b"A", 0, 8).unwrap(), Some(0));
        assert_eq!(s.find(b"CG", 0, 8).unwrap(), Some(1));
        assert_eq!(s.find(b"TTT", 0, 8).unwrap(), None);

        // empty needle => start
        assert_eq!(s.find(b"", 0, 8).unwrap(), Some(0));
        assert_eq!(s.find(b"", 3, 8).unwrap(), Some(3));

        // range-limited
        assert_eq!(s.find(b"AC", 1, 8).unwrap(), Some(4));
        assert_eq!(s.find(b"AC", 5, 8).unwrap(), None);
    }

    #[test]
    fn rfind_basic() {
        let s = DnaSeq::new(b"ACGTACGT".to_vec()).unwrap();

        assert_eq!(s.rfind(b"A", 0, 8).unwrap(), Some(4));
        assert_eq!(s.rfind(b"CG", 0, 8).unwrap(), Some(5));
        assert_eq!(s.rfind(b"TTT", 0, 8).unwrap(), None);

        // empty needle => end
        assert_eq!(s.rfind(b"", 0, 8).unwrap(), Some(8));
        assert_eq!(s.rfind(b"", 3, 8).unwrap(), Some(8));

        // range-limited
        assert_eq!(s.rfind(b"AC", 1, 8).unwrap(), Some(4));
        assert_eq!(s.rfind(b"AC", 5, 8).unwrap(), None);
        assert_eq!(s.rfind(b"AC", 0, 4).unwrap(), Some(0));
    }

    #[test]
    fn transcribe_basic() {
        let s = DnaSeq::new(b"ATGC".to_vec()).unwrap();
        let rna = s.transcribe();
        assert_eq!(rna.as_bytes(), b"AUGC");
    }

    #[test]
    fn translate_strict_rejects_non_multiple_of_3() {
        let s = DnaSeq::new(b"ATGA".to_vec()).unwrap();
        assert!(s.translate().is_err());
    }

    #[test]
    fn translate_frame_one_drops_trailing() {
        // ATGGCC + A trailing = "MA", drops 1 base
        let s = DnaSeq::new(b"ATGGCCA".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::One).unwrap();
        assert_eq!(p.as_bytes(), b"MA");
    }

    #[test]
    fn translate_frame_two() {
        // offset 1: ATGGCC -> "MA"
        let s = DnaSeq::new(b"AATGGCC".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::Two).unwrap();
        assert_eq!(p.as_bytes(), b"MA");
    }

    #[test]
    fn translate_frame_three() {
        // offset 2: ATGGCC -> "MA"
        let s = DnaSeq::new(b"CCATGGCC".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::Three).unwrap();
        assert_eq!(p.as_bytes(), b"MA");
    }

    #[test]
    fn translate_frame_short_sequences() {
        let empty = DnaSeq::new(Vec::new()).unwrap();
        assert_eq!(
            empty.translate_frame(TranslationFrame::One).unwrap().as_bytes(),
            b""
        );
        let one = DnaSeq::new(b"A".to_vec()).unwrap();
        assert_eq!(
            one.translate_frame(TranslationFrame::Two).unwrap().as_bytes(),
            b""
        );
        let two = DnaSeq::new(b"AT".to_vec()).unwrap();
        assert_eq!(
            two.translate_frame(TranslationFrame::Three).unwrap().as_bytes(),
            b""
        );
    }

    #[test]
    fn translate_frame_auto_picks_best() {
        // Frame 2 (offset 1) has ATG... -> M, frames 1 and 3 don't start with M
        // "C" + "ATGAAATTT" -> frame 2 gives "MKF"
        let s = DnaSeq::new(b"CATGAAATTT".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::Auto).unwrap();
        assert_eq!(p.as_bytes(), b"MKF");
    }

    #[test]
    fn translate_frame_auto_longest_orf() {
        // Frame 1: ATGTAA = "M*" (ORF len 1)
        // Frame 2: offset 1 -> TGT AAG CC -> "CK" (no M -> no ORF)
        // Frame 3: offset 2 -> GTA AGC C -> "VS" (no M -> no ORF)
        // Frame 1 wins (only one with M)
        let s = DnaSeq::new(b"ATGTAAGCC".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::Auto).unwrap();
        // Frame 1: ATG TAA GCC -> "M*A"
        assert_eq!(p.as_bytes(), b"M*A");
    }
}
