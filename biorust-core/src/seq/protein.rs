use crate::alphabets::protein;
use crate::error::{BioError, BioResult};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProteinSeq {
    bytes: Vec<u8>,
}

impl ProteinSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        if !protein::iupac_alphabet().is_word(bytes.as_slice()) {
            return Err(BioError::InvalidChar { ch: '?', pos: 0 });
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
}
