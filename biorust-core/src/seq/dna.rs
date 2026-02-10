use crate::alphabets::dna;
use crate::error::{BioError, BioResult};
use crate::seq::bytes::{self, IntoNeedle, Needle};
use crate::seq::protein::ProteinSeq;
use crate::seq::rna::RnaSeq;
use crate::seq::traits::SeqBytes;

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
        let mut out = Vec::with_capacity(bytes.len() / 3);
        let base_index = &*BASE_INDEX;

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

        Ok(ProteinSeq::from_bytes_unchecked(out))
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
}
