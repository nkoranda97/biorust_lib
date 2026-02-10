use crate::alphabets::rna;
use crate::error::{BioError, BioResult};
use crate::seq::bytes::{self, IntoNeedle, Needle};
use crate::seq::dna::{DnaSeq, ReverseComplement};
use crate::seq::protein::ProteinSeq;
use crate::seq::traits::SeqBytes;
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
}
