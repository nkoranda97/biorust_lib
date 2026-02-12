use crate::alphabets::rna;
use crate::error::{BioError, BioResult};
use crate::seq::bytes::{self, IntoNeedle, Needle};
use crate::seq::dna::{DnaSeq, ReverseComplement};
use crate::seq::protein::ProteinSeq;
use crate::seq::traits::SeqBytes;
use crate::seq::{best_frame_index, TranslationFrame};
use std::sync::LazyLock;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct RnaSeq {
    bytes: Vec<u8>,
}

impl RnaSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        let alphabet = rna::iupac_alphabet();
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
        let out = rna::reverse_complement(self.as_bytes());
        Self { bytes: out }
    }

    pub fn complement(&self) -> Self {
        let mut out = Vec::with_capacity(self.bytes.len());
        for &base in self.as_bytes() {
            out.push(rna::complement(base));
        }
        Self { bytes: out }
    }

    pub fn back_transcribe(&self) -> DnaSeq {
        let mut out = self.bytes.clone();
        for b in &mut out {
            if *b == b'U' {
                *b = b'T';
            } else if *b == b'u' {
                *b = b't';
            }
        }
        DnaSeq::from_bytes_unchecked(out)
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
                let idx = best_frame_index([&candidates[0], &candidates[1], &candidates[2]]);
                Ok(ProteinSeq::from_bytes_unchecked(std::mem::take(
                    &mut candidates[idx],
                )))
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

impl SeqBytes for RnaSeq {
    fn as_bytes(&self) -> &[u8] {
        RnaSeq::as_bytes(self)
    }

    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self> {
        RnaSeq::new(bytes)
    }
}

impl ReverseComplement for RnaSeq {
    fn reverse_complement(&self) -> Self {
        RnaSeq::reverse_complement(self)
    }
}

impl<'a> IntoNeedle<'a> for &'a RnaSeq {
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
    map[b'U' as usize] = 3;
    map[b'a' as usize] = 0;
    map[b'c' as usize] = 1;
    map[b'g' as usize] = 2;
    map[b'u' as usize] = 3;
    map
});

const CODON_TABLE: [u8; 64] = *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_basic() {
        let rna = RnaSeq::new(b"ACGU".to_vec()).unwrap();
        assert_eq!(rna.complement().as_bytes(), b"UGCA");
        assert_eq!(rna.reverse_complement().as_bytes(), b"ACGU");
    }

    #[test]
    fn back_transcribe_basic() {
        let rna = RnaSeq::new(b"AUGC".to_vec()).unwrap();
        let dna = rna.back_transcribe();
        assert_eq!(dna.as_bytes(), b"ATGC");
    }

    #[test]
    fn translate_basic() {
        let rna = RnaSeq::new(b"AUGGCC".to_vec()).unwrap();
        let protein = rna.translate().unwrap();
        assert_eq!(protein.as_bytes(), b"MA");
    }

    #[test]
    fn translate_strict_rejects_non_multiple_of_3() {
        let s = RnaSeq::new(b"AUGA".to_vec()).unwrap();
        assert!(s.translate().is_err());
    }

    #[test]
    fn translate_frame_one_drops_trailing() {
        let s = RnaSeq::new(b"AUGGCCA".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::One).unwrap();
        assert_eq!(p.as_bytes(), b"MA");
    }

    #[test]
    fn translate_frame_two() {
        let s = RnaSeq::new(b"AAUGGCC".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::Two).unwrap();
        assert_eq!(p.as_bytes(), b"MA");
    }

    #[test]
    fn translate_frame_auto() {
        let s = RnaSeq::new(b"CAUGGCC".to_vec()).unwrap();
        let p = s.translate_frame(TranslationFrame::Auto).unwrap();
        // Frame 2 (offset 1): AUGGCC -> "MA" has M, best ORF
        assert_eq!(p.as_bytes(), b"MA");
    }
}
