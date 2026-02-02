use crate::alphabets::dna;
use crate::error::{BioError, BioResult};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DnaSeq {
    bytes: Vec<u8>,
}

impl DnaSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        if !dna::iupac_alphabet().is_word(bytes.as_slice()) {
            return Err(BioError::InvalidChar { ch: '?', pos: 0 });
        }
        Ok(Self { bytes })
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub fn reverse_complement(&self) -> Self {
        let out = dna::reverse_complement(self.as_bytes());
        Self { bytes: out }
    }
}