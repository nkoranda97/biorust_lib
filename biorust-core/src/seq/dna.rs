use crate::alphabets::dna;
use crate::error::{BioError, BioResult};

use memchr::{memchr_iter, memmem};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DnaSeq {
    bytes: Vec<u8>,
}

impl DnaSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        if !dna::iupac_alphabet().is_word(bytes.as_slice()) {
            // keep your placeholder error for now
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

    pub fn complement(&self) -> Self {
        let mut out = Vec::with_capacity(self.bytes.len());
        for &base in self.as_bytes() {
            out.push(dna::complement(base));
        }
        Self { bytes: out }
    }

    pub fn count<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoDnaNeedle<'a>,
    {
        let hay = self.as_bytes();
        let needle = sub.into_needle()?;

        // Python: b"ABC".count(b"") == 4
        if needle.is_empty_bytes() {
            return Ok(hay.len() + 1);
        }

        match needle {
            Needle::Byte(b) => Ok(count_single_byte(hay, b)),
            Needle::Bytes(pat) => Ok(count_subslice_nonoverlapping(hay, pat)),
        }
    }

    pub fn count_overlap<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoDnaNeedle<'a>,
    {
        let hay = self.as_bytes();
        let needle = sub.into_needle()?;

        match needle {
            Needle::Byte(b) => Ok(count_single_byte(hay, b)),
            Needle::Bytes(pat) => {
                // empty matches between every char + ends => len + 1
                if pat.is_empty() {
                    return Ok(hay.len() + 1);
                }

                let finder = memmem::Finder::new(pat);

                let mut count = 0usize;
                let mut i = 0usize;

                while i <= hay.len().saturating_sub(pat.len()) {
                    if let Some(pos) = finder.find(&hay[i..]) {
                        count += 1;
                        // overlap: advance by 1 past the start of the match
                        i += pos + 1;
                    } else {
                        break;
                    }
                }

                Ok(count)
            }
        }
    }
    pub fn contains<'a, N>(&'a self, sub: N) -> BioResult<bool>
    where
        N: IntoDnaNeedle<'a>,
    {
        let hay = self.as_bytes();
        let needle = sub.into_needle()?;

        match needle {
            Needle::Byte(b) => Ok(memchr::memchr(b, hay).is_some()),

            Needle::Bytes(pat) => {
                if pat.is_empty() {
                    return Ok(true);
                }

                Ok(memmem::find(hay, pat).is_some())
            }
        }
    }

    pub fn find<'a, N>(&'a self, sub: N, start: usize, end: usize) -> BioResult<Option<usize>>
    where
        N: IntoDnaNeedle<'a>,
    {
        let hay = self.as_bytes();

        let len = hay.len();
        let start = start.min(len);
        let end = end.min(len);
        if start > end {
            return Ok(None);
        }

        let needle = sub.into_needle()?;

        match needle {
            Needle::Byte(b) => {
                let window = &hay[start..end];
                Ok(memchr::memchr(b, window).map(|i| start + i))
            }
            Needle::Bytes(pat) => {
                if pat.is_empty() {
                    return Ok(Some(start));
                }
                if pat.len() > end - start {
                    return Ok(None);
                }

                let window = &hay[start..end];
                Ok(memmem::find(window, pat).map(|i| start + i))
            }
        }
    }

    pub fn rfind<'a, N>(&'a self, sub: N, start: usize, end: usize) -> BioResult<Option<usize>>
    where
        N: IntoDnaNeedle<'a>,
    {
        let hay = self.as_bytes();

        let len = hay.len();
        let start = start.min(len);
        let end = end.min(len);
        if start > end {
            return Ok(None);
        }

        let needle = sub.into_needle()?;

        match needle {
            Needle::Byte(b) => {
                let window = &hay[start..end];
                Ok(memchr::memrchr(b, window).map(|i| start + i))
            }
            Needle::Bytes(pat) => {
                if pat.is_empty() {
                    return Ok(Some(end));
                }
                if pat.len() > end - start {
                    return Ok(None);
                }

                let window = &hay[start..end];
                let finder = memmem::Finder::new(pat);

                let mut pos = 0usize;
                let mut last = None;

                while pos <= window.len().saturating_sub(pat.len()) {
                    match finder.find(&window[pos..]) {
                        Some(i) => {
                            let found = pos + i;
                            last = Some(found);
                            pos = found + 1;
                        }
                        None => break,
                    }
                }

                Ok(last.map(|i| start + i))
            }
        }
    }
}

/// Internal “needle” representation.
#[derive(Copy, Clone, Debug)]
pub enum Needle<'a> {
    Bytes(&'a [u8]),
    Byte(u8),
}

impl Needle<'_> {
    #[inline]
    fn is_empty_bytes(&self) -> bool {
        matches!(self, Needle::Bytes(b) if b.is_empty())
    }
}

pub trait IntoDnaNeedle<'a> {
    fn into_needle(self) -> BioResult<Needle<'a>>;
}

impl<'a> IntoDnaNeedle<'a> for &'a [u8] {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self))
    }
}

impl<'a> IntoDnaNeedle<'a> for &'a str {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_bytes()))
    }
}

impl<'a> IntoDnaNeedle<'a> for u8 {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(self))
    }
}

impl<'a> IntoDnaNeedle<'a> for &'a DnaSeq {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_bytes()))
    }
}

impl<'a, const N: usize> IntoDnaNeedle<'a> for &'a [u8; N] {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_slice()))
    }
}

#[inline]
fn checked_int_byte(x: i128) -> BioResult<u8> {
    if (0..=255).contains(&x) {
        Ok(x as u8)
    } else {
        Err(BioError::IntByteOutOfRange { val: x })
    }
}

impl<'a> IntoDnaNeedle<'a> for i64 {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(checked_int_byte(self as i128)?))
    }
}

impl<'a> IntoDnaNeedle<'a> for isize {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(checked_int_byte(self as i128)?))
    }
}

impl<'a> IntoDnaNeedle<'a> for usize {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(checked_int_byte(self as i128)?))
    }
}

#[inline]
fn count_single_byte(hay: &[u8], b: u8) -> usize {
    memchr_iter(b, hay).count()
}

fn count_subslice_nonoverlapping(hay: &[u8], needle: &[u8]) -> usize {
    debug_assert!(!needle.is_empty());

    if needle.len() == 1 {
        return count_single_byte(hay, needle[0]);
    }

    let finder = memmem::Finder::new(needle);
    let mut count = 0usize;
    let mut pos = 0usize;

    while pos <= hay.len() {
        match finder.find(&hay[pos..]) {
            Some(i) => {
                count += 1;
                pos += i + needle.len(); // non-overlapping
            }
            None => break,
        }
    }

    count
}

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
}
