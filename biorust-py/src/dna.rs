#![allow(clippy::useless_conversion)]

use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PySlice, PyString, PyTuple};

use crate::protein::Protein;
use crate::utils::{self, PyDnaNeedle};
use biorust_core::seq::dna::DnaSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
pub struct DNA {
    pub(crate) inner: DnaSeq,
}

#[pymethods]
impl DNA {
    #[new]
    fn new(seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let bytes: Vec<u8> = if let Ok(s) = seq.downcast::<PyString>() {
            s.to_str()?.as_bytes().to_vec()
        } else {
            seq.extract::<Vec<u8>>()
                .map_err(|_| PyValueError::new_err("DNA() expects str or bytes-like input"))?
        };

        let inner = DnaSeq::new(bytes).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn reverse_complement(&self) -> Self {
        Self {
            inner: self.inner.reverse_complement(),
        }
    }

    fn complement(&self) -> Self {
        Self {
            inner: self.inner.complement(),
        }
    }

    fn translate(&self) -> Protein {
        Protein {
            inner: self.inner.translate(),
        }
    }

    #[inline]
    pub(crate) fn as_bytes(&self) -> &[u8] {
        self.inner.as_bytes()
    }

    fn to_bytes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, self.as_bytes())
    }

    fn __len__(&self) -> usize {
        self.as_bytes().len()
    }

    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        let other = utils::extract_bytes(other)
            .map_err(|_| PyTypeError::new_err("expected DNA, str, or bytes-like object"))?;

        match op {
            CompareOp::Eq => Ok(self.as_bytes() == other.as_slice()),
            CompareOp::Ne => Ok(self.as_bytes() != other.as_slice()),
            CompareOp::Lt => Ok(self.as_bytes() < other.as_slice()),
            CompareOp::Le => Ok(self.as_bytes() <= other.as_slice()),
            CompareOp::Gt => Ok(self.as_bytes() > other.as_slice()),
            CompareOp::Ge => Ok(self.as_bytes() >= other.as_slice()),
        }
    }

    fn __bytes__<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, self.as_bytes())
    }

    fn __str__(&self) -> PyResult<String> {
        std::str::from_utf8(self.as_bytes())
            .map(|s| s.to_string())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> PyResult<String> {
        let s = std::str::from_utf8(self.inner.as_bytes()).unwrap_or("<bytes>");
        Ok(format!("DNA({s:?})"))
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let bytes = self.inner.as_bytes();

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
                    i += step; // negative
                }
            }

            let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            return Ok(Py::new(py, DNA { inner })?.to_object(py));
        }

        let index: isize = index
            .extract()
            .map_err(|_| PyTypeError::new_err("index must be int or slice"))?;

        let n = bytes.len() as isize;
        let i = if index < 0 { index + n } else { index };

        if i < 0 || i >= n {
            return Err(PyIndexError::new_err("index out of range"));
        }

        let one = vec![bytes[i as usize]];
        let inner = DnaSeq::new(one).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Py::new(py, DNA { inner })?.to_object(py))
    }

    fn __radd__(&self, seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let mut left = utils::extract_bytes(seq)?;

        // left + self
        left.extend_from_slice(self.as_bytes());

        let inner = DnaSeq::new(left).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn __add__(&self, seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let right = utils::extract_bytes(seq)?;

        // self + right
        let mut out = Vec::with_capacity(self.inner.as_bytes().len() + right.len());
        out.extend_from_slice(self.as_bytes());
        out.extend_from_slice(&right);

        let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn __mul__(&self, num: &Bound<'_, PyAny>) -> PyResult<Self> {
        let n: isize = num
            .extract()
            .map_err(|_| PyTypeError::new_err("num must be int"))?;

        if n <= 0 {
            let inner =
                DnaSeq::new(Vec::new()).map_err(|e| PyValueError::new_err(e.to_string()))?;
            return Ok(Self { inner });
        }

        let n: usize = n as usize;
        let bytes = self.inner.as_bytes();

        let total_len = bytes
            .len()
            .checked_mul(n)
            .ok_or_else(|| PyValueError::new_err("resulting sequence is too large"))?;

        let mut out = Vec::with_capacity(total_len);
        for _ in 0..n {
            out.extend_from_slice(bytes);
        }

        let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn __rmul__(&self, num: &Bound<'_, PyAny>) -> PyResult<Self> {
        self.__mul__(num)
    }

    fn count(&self, sub: &Bound<'_, PyAny>) -> PyResult<usize> {
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.count(&other.inner),
            PyDnaNeedle::Bytes(bytes) => self.inner.count(bytes.as_slice()),
            PyDnaNeedle::Byte(b) => self.inner.count(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn count_overlap(&self, sub: &Bound<'_, PyAny>) -> PyResult<usize> {
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.count_overlap(&other.inner),
            PyDnaNeedle::Bytes(bytes) => self.inner.count_overlap(bytes.as_slice()),
            PyDnaNeedle::Byte(b) => self.inner.count_overlap(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __contains__(&self, sub: &Bound<'_, PyAny>) -> PyResult<bool> {
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.contains(&other.inner),
            PyDnaNeedle::Bytes(bytes) => self.inner.contains(bytes.as_slice()),
            PyDnaNeedle::Byte(b) => self.inner.contains(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[pyo3(signature = (prefix, start=None, end=None))]
    fn startswith(
        &self,
        prefix: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<bool> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let bytes = self.as_bytes();
        let window: &[u8] = if s <= e { &bytes[s..e] } else { &bytes[0..0] };

        let matches = |needle: PyDnaNeedle<'_>| -> bool {
            match needle {
                PyDnaNeedle::Dna(other) => window.starts_with(other.as_bytes()),
                PyDnaNeedle::Bytes(seq) => window.starts_with(seq.as_slice()),
                PyDnaNeedle::Byte(b) => window.first().map(|v| *v == b).unwrap_or(false),
            }
        };

        if let Ok(tuple) = prefix.downcast::<PyTuple>() {
            for item in tuple.iter() {
                let needle = utils::extract_dna_needle(&item)?;
                if matches(needle) {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        let needle = utils::extract_dna_needle(prefix)?;
        Ok(matches(needle))
    }

    #[pyo3(signature = (suffix, start=None, end=None))]
    fn endswith(
        &self,
        suffix: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<bool> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let bytes = self.as_bytes();
        let window: &[u8] = if s <= e { &bytes[s..e] } else { &bytes[0..0] };

        let matches = |needle: PyDnaNeedle<'_>| -> bool {
            match needle {
                PyDnaNeedle::Dna(other) => window.ends_with(other.as_bytes()),
                PyDnaNeedle::Bytes(seq) => window.ends_with(seq.as_slice()),
                PyDnaNeedle::Byte(b) => window.last().map(|v| *v == b).unwrap_or(false),
            }
        };

        if let Ok(tuple) = suffix.downcast::<PyTuple>() {
            for item in tuple.iter() {
                let needle = utils::extract_dna_needle(&item)?;
                if matches(needle) {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        let needle = utils::extract_dna_needle(suffix)?;
        Ok(matches(needle))
    }

    #[pyo3(signature = (sep=None, maxsplit=-1))]
    fn split<'py>(
        &self,
        py: Python<'py>,
        sep: Option<&Bound<'py, PyAny>>,
        maxsplit: isize,
    ) -> PyResult<Vec<Py<DNA>>> {
        let bytes = self.as_bytes();
        let parts = match sep {
            None => split_on_whitespace(bytes, maxsplit),
            Some(obj) => split_on_sep(bytes, obj, maxsplit)?,
        };

        dna_list_from_parts(py, parts)
    }

    #[pyo3(signature = (sep=None, maxsplit=-1))]
    fn rsplit<'py>(
        &self,
        py: Python<'py>,
        sep: Option<&Bound<'py, PyAny>>,
        maxsplit: isize,
    ) -> PyResult<Vec<Py<DNA>>> {
        let bytes = self.as_bytes();
        let parts = match sep {
            None => rsplit_on_whitespace(bytes, maxsplit),
            Some(obj) => rsplit_on_sep(bytes, obj, maxsplit)?,
        };

        dna_list_from_parts(py, parts)
    }

    #[pyo3(signature = (chars=None))]
    fn strip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let (start, end) = trim_range(bytes, chars, true, true)?;
        let inner = DnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn lstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let (start, end) = trim_range(bytes, chars, true, false)?;
        let inner = DnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn rstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let (start, end) = trim_range(bytes, chars, false, true)?;
        let inner = DnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn upper(&self) -> PyResult<Self> {
        let out: Vec<u8> = self
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn lower(&self) -> PyResult<Self> {
        let out: Vec<u8> = self
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_lowercase())
            .collect();
        let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (sub, start=None, end=None))]
    fn find(
        &self,
        sub: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<isize> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.find(&other.inner, s, e),
            PyDnaNeedle::Bytes(bytes) => self.inner.find(bytes.as_slice(), s, e),
            PyDnaNeedle::Byte(b) => self.inner.find(b, s, e),
        };

        match res.map_err(|err| PyValueError::new_err(err.to_string()))? {
            Some(pos) => Ok(pos as isize),
            None => Ok(-1),
        }
    }

    #[pyo3(signature = (sub, start=None, end=None))]
    fn index(
        &self,
        sub: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<isize> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.find(&other.inner, s, e),
            PyDnaNeedle::Bytes(bytes) => self.inner.find(bytes.as_slice(), s, e),
            PyDnaNeedle::Byte(b) => self.inner.find(b, s, e),
        };

        match res.map_err(|err| PyValueError::new_err(err.to_string()))? {
            Some(pos) => Ok(pos as isize),
            None => Err(PyValueError::new_err("subsection not found")),
        }
    }

    #[pyo3(signature = (sub, start=None, end=None))]
    fn rfind(
        &self,
        sub: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<isize> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.rfind(&other.inner, s, e),
            PyDnaNeedle::Bytes(bytes) => self.inner.rfind(bytes.as_slice(), s, e),
            PyDnaNeedle::Byte(b) => self.inner.rfind(b, s, e),
        };

        match res.map_err(|err| PyValueError::new_err(err.to_string()))? {
            Some(pos) => Ok(pos as isize),
            None => Ok(-1),
        }
    }

    #[pyo3(signature = (sub, start=None, end=None))]
    fn rindex(
        &self,
        sub: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<isize> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let needle = utils::extract_dna_needle(sub)?;

        let res = match needle {
            PyDnaNeedle::Dna(other) => self.inner.rfind(&other.inner, s, e),
            PyDnaNeedle::Bytes(bytes) => self.inner.rfind(bytes.as_slice(), s, e),
            PyDnaNeedle::Byte(b) => self.inner.rfind(b, s, e),
        };

        match res.map_err(|err| PyValueError::new_err(err.to_string()))? {
            Some(pos) => Ok(pos as isize),
            None => Err(PyValueError::new_err("subsection not found")),
        }
    }
}

#[pyfunction]
fn complement(seq: &Bound<'_, PyAny>) -> PyResult<DNA> {
    if let Ok(dna) = seq.extract::<PyRef<'_, DNA>>() {
        return Ok(DNA {
            inner: dna.inner.complement(),
        });
    }

    let bytes = utils::extract_bytes(seq)?;
    let inner = DnaSeq::new(bytes).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(DNA {
        inner: inner.complement(),
    })
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DNA>()?;
    m.add_function(wrap_pyfunction!(complement, m)?)?;
    Ok(())
}

fn dna_list_from_parts(py: Python<'_>, parts: Vec<Vec<u8>>) -> PyResult<Vec<Py<DNA>>> {
    let mut out = Vec::with_capacity(parts.len());
    for part in parts {
        let inner = DnaSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
        out.push(Py::new(py, DNA { inner })?);
    }
    Ok(out)
}

fn split_on_whitespace(hay: &[u8], maxsplit: isize) -> Vec<Vec<u8>> {
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

fn rsplit_on_whitespace(hay: &[u8], maxsplit: isize) -> Vec<Vec<u8>> {
    let len = hay.len();
    let maxsplit = if maxsplit < 0 {
        usize::MAX
    } else {
        maxsplit as usize
    };

    let mut out = Vec::new();
    let mut i = len;

    if maxsplit == 0 {
        while i > 0 && hay[i - 1].is_ascii_whitespace() {
            i -= 1;
        }
        if i == 0 {
            return out;
        }
        let mut start = 0usize;
        while start < i && hay[start].is_ascii_whitespace() {
            start += 1;
        }
        out.push(hay[start..i].to_vec());
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

fn split_on_sep(hay: &[u8], sep: &Bound<'_, PyAny>, maxsplit: isize) -> PyResult<Vec<Vec<u8>>> {
    let needle = utils::extract_dna_needle(sep)?;

    match needle {
        PyDnaNeedle::Byte(b) => Ok(split_on_byte(hay, b, maxsplit)),
        PyDnaNeedle::Bytes(bytes) => split_on_bytes(hay, bytes.as_slice(), maxsplit),
        PyDnaNeedle::Dna(other) => split_on_bytes(hay, other.as_bytes(), maxsplit),
    }
}

fn rsplit_on_sep(hay: &[u8], sep: &Bound<'_, PyAny>, maxsplit: isize) -> PyResult<Vec<Vec<u8>>> {
    let needle = utils::extract_dna_needle(sep)?;

    match needle {
        PyDnaNeedle::Byte(b) => Ok(rsplit_on_byte(hay, b, maxsplit)),
        PyDnaNeedle::Bytes(bytes) => rsplit_on_bytes(hay, bytes.as_slice(), maxsplit),
        PyDnaNeedle::Dna(other) => rsplit_on_bytes(hay, other.as_bytes(), maxsplit),
    }
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

fn trim_range(
    hay: &[u8],
    chars: Option<&Bound<'_, PyAny>>,
    left: bool,
    right: bool,
) -> PyResult<(usize, usize)> {
    let len = hay.len();
    let mut start = 0usize;
    let mut end = len;

    let mut mask = [false; 256];
    let mut use_mask = false;
    let mut single_byte: Option<u8> = None;

    if let Some(obj) = chars {
        let needle = utils::extract_dna_needle(obj)?;
        let bytes = match needle {
            PyDnaNeedle::Dna(other) => other.as_bytes().to_vec(),
            PyDnaNeedle::Bytes(bytes) => bytes,
            PyDnaNeedle::Byte(b) => vec![b],
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
