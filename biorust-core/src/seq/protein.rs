use crate::alphabets::protein;
use crate::error::{BioError, BioResult};
use crate::seq::bytes::{self, IntoNeedle, Needle};

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

    pub fn count<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::count(self.as_bytes(), needle))
    }

    pub fn count_overlap<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::count_overlap(self.as_bytes(), needle))
    }

    pub fn contains<'a, N>(&'a self, sub: N) -> BioResult<bool>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::contains(self.as_bytes(), needle))
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

impl<'a> IntoNeedle<'a> for &'a ProteinSeq {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_bytes()))
    }
}
