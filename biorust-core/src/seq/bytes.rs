use crate::error::{BioError, BioResult};

use memchr::{memchr_iter, memmem};

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

pub trait IntoNeedle<'a> {
    fn into_needle(self) -> BioResult<Needle<'a>>;
}

impl<'a> IntoNeedle<'a> for &'a [u8] {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self))
    }
}

impl<'a> IntoNeedle<'a> for &'a str {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_bytes()))
    }
}

impl<'a> IntoNeedle<'a> for u8 {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(self))
    }
}

impl<'a, const N: usize> IntoNeedle<'a> for &'a [u8; N] {
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

impl<'a> IntoNeedle<'a> for i64 {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(checked_int_byte(self as i128)?))
    }
}

impl<'a> IntoNeedle<'a> for isize {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(checked_int_byte(self as i128)?))
    }
}

impl<'a> IntoNeedle<'a> for usize {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Byte(checked_int_byte(self as i128)?))
    }
}

pub fn count(hay: &[u8], needle: Needle<'_>) -> usize {
    // Python: b"ABC".count(b"") == 4
    if needle.is_empty_bytes() {
        return hay.len() + 1;
    }

    match needle {
        Needle::Byte(b) => count_single_byte(hay, b),
        Needle::Bytes(pat) => count_subslice_nonoverlapping(hay, pat),
    }
}

pub fn count_overlap(hay: &[u8], needle: Needle<'_>) -> usize {
    match needle {
        Needle::Byte(b) => count_single_byte(hay, b),
        Needle::Bytes(pat) => {
            if pat.is_empty() {
                return hay.len() + 1;
            }

            let finder = memmem::Finder::new(pat);
            let mut count = 0usize;
            let mut i = 0usize;

            while i <= hay.len().saturating_sub(pat.len()) {
                if let Some(pos) = finder.find(&hay[i..]) {
                    count += 1;
                    i += pos + 1;
                } else {
                    break;
                }
            }

            count
        }
    }
}

pub fn contains(hay: &[u8], needle: Needle<'_>) -> bool {
    match needle {
        Needle::Byte(b) => memchr::memchr(b, hay).is_some(),
        Needle::Bytes(pat) => {
            if pat.is_empty() {
                return true;
            }
            memmem::find(hay, pat).is_some()
        }
    }
}

pub fn find(hay: &[u8], needle: Needle<'_>, start: usize, end: usize) -> Option<usize> {
    let len = hay.len();
    let start = start.min(len);
    let end = end.min(len);
    if start > end {
        return None;
    }

    match needle {
        Needle::Byte(b) => {
            let window = &hay[start..end];
            memchr::memchr(b, window).map(|i| start + i)
        }
        Needle::Bytes(pat) => {
            if pat.is_empty() {
                return Some(start);
            }
            if pat.len() > end - start {
                return None;
            }
            let window = &hay[start..end];
            memmem::find(window, pat).map(|i| start + i)
        }
    }
}

pub fn rfind(hay: &[u8], needle: Needle<'_>, start: usize, end: usize) -> Option<usize> {
    let len = hay.len();
    let start = start.min(len);
    let end = end.min(len);
    if start > end {
        return None;
    }

    match needle {
        Needle::Byte(b) => {
            let window = &hay[start..end];
            memchr::memrchr(b, window).map(|i| start + i)
        }
        Needle::Bytes(pat) => {
            if pat.is_empty() {
                return Some(end);
            }
            if pat.len() > end - start {
                return None;
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

            last.map(|i| start + i)
        }
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
                pos += i + needle.len();
            }
            None => break,
        }
    }

    count
}
