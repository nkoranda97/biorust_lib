use crate::alphabets::dna;
use crate::error::{BioError, BioResult};
use crate::seq::dna::DnaSeq;
use crate::seq::traits::SeqBytes;

use std::sync::LazyLock;

static GAPPED_DNA_IUPAC: LazyLock<bit_set::BitSet> = LazyLock::new(|| {
    let mut s = dna::iupac_alphabet().symbols;
    s.insert(b'-' as usize);
    s.insert(b'.' as usize);
    s
});

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct GappedDnaSeq {
    bytes: Vec<u8>,
}

impl GappedDnaSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        let symbols = &*GAPPED_DNA_IUPAC;
        for (pos, &b) in bytes.iter().enumerate() {
            if !symbols.contains(b as usize) {
                return Err(BioError::InvalidChar { ch: b as char, pos });
            }
        }
        Ok(Self { bytes })
    }

    #[inline]
    #[allow(dead_code)]
    pub(crate) fn from_bytes_unchecked(bytes: Vec<u8>) -> Self {
        Self { bytes }
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Strip gap characters (`-` and `.`) and return a strict `DnaSeq`.
    pub fn ungapped(&self) -> DnaSeq {
        let out: Vec<u8> = self
            .bytes
            .iter()
            .copied()
            .filter(|&b| b != b'-' && b != b'.')
            .collect();
        DnaSeq::from_bytes_unchecked(out)
    }
}

impl SeqBytes for GappedDnaSeq {
    fn as_bytes(&self) -> &[u8] {
        GappedDnaSeq::as_bytes(self)
    }

    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self> {
        GappedDnaSeq::new(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_gapped_seq() {
        let seq = GappedDnaSeq::new(b"AC-GT.N".to_vec()).unwrap();
        assert_eq!(seq.as_bytes(), b"AC-GT.N");
        assert_eq!(seq.len(), 7);
        assert!(!seq.is_empty());
    }

    #[test]
    fn invalid_char_rejected() {
        let err = GappedDnaSeq::new(b"AC#GT".to_vec()).unwrap_err();
        match err {
            BioError::InvalidChar { ch, pos } => {
                assert_eq!(ch, '#');
                assert_eq!(pos, 2);
            }
            _ => panic!("expected InvalidChar"),
        }
    }

    #[test]
    fn ungap_round_trip() {
        let gapped = GappedDnaSeq::new(b"A-C.G-T".to_vec()).unwrap();
        let ungapped = gapped.ungapped();
        assert_eq!(ungapped.as_bytes(), b"ACGT");
    }

    #[test]
    fn empty_seq() {
        let seq = GappedDnaSeq::new(Vec::new()).unwrap();
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
        assert_eq!(seq.ungapped().len(), 0);
    }

    #[test]
    fn all_gaps() {
        let seq = GappedDnaSeq::new(b"---...".to_vec()).unwrap();
        assert_eq!(seq.len(), 6);
        assert_eq!(seq.ungapped().len(), 0);
    }

    #[test]
    fn iupac_chars_accepted() {
        // All IUPAC + gaps should pass
        let seq = GappedDnaSeq::new(b"ACGTRYSWKMBDHVNacgtryswkmbdhvn-.".to_vec()).unwrap();
        assert_eq!(seq.len(), 32);
    }
}
