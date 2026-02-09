use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyBytes, PySlice};
use pyo3::PyClass;

use crate::utils::normalize_range;

pub fn seq_to_bytes<'py>(py: Python<'py>, bytes: &[u8]) -> Bound<'py, PyBytes> {
    PyBytes::new_bound(py, bytes)
}

pub fn seq_str(bytes: &[u8]) -> PyResult<String> {
    std::str::from_utf8(bytes)
        .map(|s| s.to_string())
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

pub fn seq_repr(bytes: &[u8], name: &str) -> String {
    let s = std::str::from_utf8(bytes).unwrap_or("<bytes>");
    format!("{name}({s:?})")
}

pub fn seq_getitem<F>(bytes: &[u8], index: &Bound<'_, PyAny>, make: F) -> PyResult<PyObject>
where
    F: Fn(Vec<u8>) -> PyResult<PyObject>,
{
    if let Ok(slice) = index.downcast::<PySlice>() {
        let idx = slice.indices(bytes.len() as isize)?;
        let (start, stop, step) = (idx.start, idx.stop, idx.step);

        let mut out = Vec::new();

        if step > 0 {
            let mut i = start;
            while i < stop {
                out.push(bytes[i as usize]);
                i += step;
            }
        } else {
            let mut i = start;
            while i > stop {
                out.push(bytes[i as usize]);
                i += step;
            }
        }

        return make(out);
    }

    let index: isize = index
        .extract()
        .map_err(|_| PyTypeError::new_err("index must be int or slice"))?;

    let n = bytes.len() as isize;
    let i = if index < 0 { index + n } else { index };

    if i < 0 || i >= n {
        return Err(PyIndexError::new_err("index out of range"));
    }

    make(vec![bytes[i as usize]])
}

pub enum NeedleBytes<'a> {
    Bytes(&'a [u8]),
    Byte(u8),
}

pub fn seq_upper<T, F>(bytes: &[u8], make: F) -> PyResult<T>
where
    F: Fn(Vec<u8>) -> PyResult<T>,
{
    let out: Vec<u8> = bytes.iter().map(|b| b.to_ascii_uppercase()).collect();
    make(out)
}

pub fn seq_lower<T, F>(bytes: &[u8], make: F) -> PyResult<T>
where
    F: Fn(Vec<u8>) -> PyResult<T>,
{
    let out: Vec<u8> = bytes.iter().map(|b| b.to_ascii_lowercase()).collect();
    make(out)
}

pub fn needle_starts_with(window: &[u8], needle: NeedleBytes<'_>) -> bool {
    match needle {
        NeedleBytes::Bytes(seq) => window.starts_with(seq),
        NeedleBytes::Byte(b) => window.first().map(|v| *v == b).unwrap_or(false),
    }
}

pub fn needle_ends_with(window: &[u8], needle: NeedleBytes<'_>) -> bool {
    match needle {
        NeedleBytes::Bytes(seq) => window.ends_with(seq),
        NeedleBytes::Byte(b) => window.last().map(|v| *v == b).unwrap_or(false),
    }
}

pub fn split_on_whitespace(hay: &[u8], maxsplit: isize) -> Vec<Vec<u8>> {
    let len = hay.len();
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    let mut out = Vec::new();
    let mut i = 0usize;

    if maxsplit == 0 {
        while i < len && hay[i].is_ascii_whitespace() {
            i += 1;
        }
        if i < len {
            out.push(hay[i..].to_vec());
        }
        return out;
    }

    let mut splits = 0usize;
    while i < len {
        while i < len && hay[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len {
            break;
        }

        if splits == maxsplit {
            out.push(hay[i..].to_vec());
            return out;
        }

        let start = i;
        while i < len && !hay[i].is_ascii_whitespace() {
            i += 1;
        }
        out.push(hay[start..i].to_vec());
        splits += 1;
    }

    out
}

pub fn rsplit_on_whitespace(hay: &[u8], maxsplit: isize) -> Vec<Vec<u8>> {
    let len = hay.len();
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    let mut out = Vec::new();
    let mut i = len;

    if maxsplit == 0 {
        // Python rsplit(maxsplit=0) strips trailing whitespace only
        while i > 0 && hay[i - 1].is_ascii_whitespace() {
            i -= 1;
        }
        if i > 0 {
            out.push(hay[..i].to_vec());
        }
        return out;
    }

    let mut splits = 0usize;
    while i > 0 {
        while i > 0 && hay[i - 1].is_ascii_whitespace() {
            i -= 1;
        }
        if i == 0 {
            break;
        }

        if splits == maxsplit {
            let mut start = 0usize;
            while start < i && hay[start].is_ascii_whitespace() {
                start += 1;
            }
            out.push(hay[start..i].to_vec());
            out.reverse();
            return out;
        }

        let end = i;
        while i > 0 && !hay[i - 1].is_ascii_whitespace() {
            i -= 1;
        }
        out.push(hay[i..end].to_vec());
        splits += 1;
    }

    out.reverse();
    out
}

pub fn split_on_sep(
    hay: &[u8],
    needle: NeedleBytes<'_>,
    maxsplit: isize,
) -> PyResult<Vec<Vec<u8>>> {
    match needle {
        NeedleBytes::Byte(b) => Ok(split_on_byte(hay, b, maxsplit)),
        NeedleBytes::Bytes(bytes) => split_on_bytes(hay, bytes, maxsplit),
    }
}

pub fn rsplit_on_sep(
    hay: &[u8],
    needle: NeedleBytes<'_>,
    maxsplit: isize,
) -> PyResult<Vec<Vec<u8>>> {
    match needle {
        NeedleBytes::Byte(b) => Ok(rsplit_on_byte(hay, b, maxsplit)),
        NeedleBytes::Bytes(bytes) => rsplit_on_bytes(hay, bytes, maxsplit),
    }
}

pub fn trim_range(
    hay: &[u8],
    chars: Option<NeedleBytes<'_>>,
    left: bool,
    right: bool,
) -> PyResult<(usize, usize)> {
    let len = hay.len();
    let mut start = 0usize;
    let mut end = len;

    let mut mask = [false; 256];
    let mut use_mask = false;
    let mut single_byte: Option<u8> = None;

    if let Some(needle) = chars {
        let bytes = match needle {
            NeedleBytes::Bytes(bytes) => bytes.to_vec(),
            NeedleBytes::Byte(b) => vec![b],
        };

        if bytes.is_empty() {
            return Ok((0, len));
        }

        if bytes.len() == 1 {
            single_byte = Some(bytes[0]);
        } else {
            use_mask = true;
            for &b in bytes.iter() {
                mask[b as usize] = true;
            }
        }
    }

    let is_trim = |b: u8, single_byte: Option<u8>, use_mask: bool, mask: &[bool; 256]| -> bool {
        if let Some(sb) = single_byte {
            return b == sb;
        }
        if use_mask {
            return mask[b as usize];
        }
        b.is_ascii_whitespace()
    };

    if left {
        while start < end && is_trim(hay[start], single_byte, use_mask, &mask) {
            start += 1;
        }
    }

    if right {
        while end > start && is_trim(hay[end - 1], single_byte, use_mask, &mask) {
            end -= 1;
        }
    }

    Ok((start, end))
}

pub fn list_from_parts<T, F>(parts: Vec<Vec<u8>>, make: F) -> PyResult<Vec<Py<T>>>
where
    T: PyClass,
    F: Fn(Vec<u8>) -> PyResult<Py<T>>,
{
    let mut out = Vec::with_capacity(parts.len());
    for part in parts {
        out.push(make(part)?);
    }
    Ok(out)
}

fn split_on_bytes(hay: &[u8], sep: &[u8], maxsplit: isize) -> PyResult<Vec<Vec<u8>>> {
    if sep.is_empty() {
        return Err(PyValueError::new_err("empty separator"));
    }

    if sep.len() == 1 {
        return Ok(split_on_byte(hay, sep[0], maxsplit));
    }

    Ok(split_on_bytes_multi(hay, sep, maxsplit))
}

fn rsplit_on_bytes(hay: &[u8], sep: &[u8], maxsplit: isize) -> PyResult<Vec<Vec<u8>>> {
    if sep.is_empty() {
        return Err(PyValueError::new_err("empty separator"));
    }

    if sep.len() == 1 {
        return Ok(rsplit_on_byte(hay, sep[0], maxsplit));
    }

    Ok(rsplit_on_bytes_multi(hay, sep, maxsplit))
}

fn split_on_byte(hay: &[u8], b: u8, maxsplit: isize) -> Vec<Vec<u8>> {
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    if maxsplit == 0 {
        return vec![hay.to_vec()];
    }

    let mut out = Vec::new();
    let mut start = 0usize;
    let mut splits = 0usize;

    for (i, &c) in hay.iter().enumerate() {
        if c == b && splits < maxsplit {
            out.push(hay[start..i].to_vec());
            start = i + 1;
            splits += 1;
        }
    }

    out.push(hay[start..].to_vec());
    out
}

fn rsplit_on_byte(hay: &[u8], b: u8, maxsplit: isize) -> Vec<Vec<u8>> {
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    if maxsplit == 0 {
        return vec![hay.to_vec()];
    }

    let mut out = Vec::new();
    let mut end = hay.len();
    let mut splits = 0usize;

    let mut i = hay.len();
    while i > 0 {
        i -= 1;
        if hay[i] == b && splits < maxsplit {
            out.push(hay[i + 1..end].to_vec());
            end = i;
            splits += 1;
        }
    }

    out.push(hay[..end].to_vec());
    out.reverse();
    out
}

fn split_on_bytes_multi(hay: &[u8], sep: &[u8], maxsplit: isize) -> Vec<Vec<u8>> {
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    if maxsplit == 0 {
        return vec![hay.to_vec()];
    }

    let mut out = Vec::new();
    let mut start = 0usize;
    let mut splits = 0usize;

    while splits < maxsplit {
        let pos = find_subslice(hay, sep, start);
        match pos {
            Some(i) => {
                out.push(hay[start..i].to_vec());
                start = i + sep.len();
                splits += 1;
            }
            None => break,
        }
    }

    out.push(hay[start..].to_vec());
    out
}

fn rsplit_on_bytes_multi(hay: &[u8], sep: &[u8], maxsplit: isize) -> Vec<Vec<u8>> {
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    if maxsplit == 0 {
        return vec![hay.to_vec()];
    }

    let mut out = Vec::new();
    let mut end = hay.len();
    let mut splits = 0usize;

    while splits < maxsplit {
        let pos = rfind_subslice(hay, sep, end);
        match pos {
            Some(i) => {
                out.push(hay[i + sep.len()..end].to_vec());
                end = i;
                splits += 1;
            }
            None => break,
        }
    }

    out.push(hay[..end].to_vec());
    out.reverse();
    out
}

fn find_subslice(hay: &[u8], needle: &[u8], start: usize) -> Option<usize> {
    if needle.len() > hay.len().saturating_sub(start) {
        return None;
    }
    hay[start..]
        .windows(needle.len())
        .position(|w| w == needle)
        .map(|i| start + i)
}

fn rfind_subslice(hay: &[u8], needle: &[u8], end: usize) -> Option<usize> {
    if needle.len() > end {
        return None;
    }
    hay[..end].windows(needle.len()).rposition(|w| w == needle)
}

pub fn startswith_window(bytes: &[u8], start: Option<isize>, end: Option<isize>) -> &[u8] {
    let (s, e) = normalize_range(bytes.len(), start, end);
    if s <= e {
        &bytes[s..e]
    } else {
        &bytes[0..0]
    }
}
