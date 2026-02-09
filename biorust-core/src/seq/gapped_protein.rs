use crate::alphabets::protein;
use crate::error::{BioError, BioResult};
use crate::seq::protein::ProteinSeq;
use crate::seq::traits::SeqBytes;

use std::sync::LazyLock;

static GAPPED_PROTEIN_IUPAC: LazyLock<bit_set::BitSet> = LazyLock::new(|| {
    let mut s = protein::iupac_alphabet().symbols;
    s.insert(b'-' as usize);
    s.insert(b'.' as usize);
    s
});

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct GappedProteinSeq {
    bytes: Vec<u8>,
}

impl GappedProteinSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        let symbols = &*GAPPED_PROTEIN_IUPAC;
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

    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Strip gap characters (`-` and `.`) and return a strict `ProteinSeq`.
    pub fn ungapped(&self) -> ProteinSeq {
        let out: Vec<u8> = self
            .bytes
            .iter()
            .copied()
            .filter(|&b| b != b'-' && b != b'.')
            .collect();
        ProteinSeq::from_bytes_unchecked(out)
    }
}

impl SeqBytes for GappedProteinSeq {
    fn as_bytes(&self) -> &[u8] {
        GappedProteinSeq::as_bytes(self)
    }

    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self> {
        GappedProteinSeq::new(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_gapped_seq() {
        let seq = GappedProteinSeq::new(b"AC-DE.F".to_vec()).unwrap();
        assert_eq!(seq.as_bytes(), b"AC-DE.F");
        assert_eq!(seq.len(), 7);
        assert!(!seq.is_empty());
    }

    #[test]
    fn invalid_char_rejected() {
        let err = GappedProteinSeq::new(b"AC#DE".to_vec()).unwrap_err();
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
        let gapped = GappedProteinSeq::new(b"A-C.D-E".to_vec()).unwrap();
        let ungapped = gapped.ungapped();
        assert_eq!(ungapped.as_bytes(), b"ACDE");
    }

    #[test]
    fn empty_seq() {
        let seq = GappedProteinSeq::new(Vec::new()).unwrap();
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
        assert_eq!(seq.ungapped().len(), 0);
    }

    #[test]
    fn all_gaps() {
        let seq = GappedProteinSeq::new(b"---...".to_vec()).unwrap();
        assert_eq!(seq.len(), 6);
        assert_eq!(seq.ungapped().len(), 0);
    }

    #[test]
    fn iupac_chars_accepted() {
        let seq =
            GappedProteinSeq::new(b"ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz*-.".to_vec())
                .unwrap();
        assert_eq!(seq.len(), 49);
    }
}
