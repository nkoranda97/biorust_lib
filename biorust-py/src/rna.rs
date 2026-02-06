#![allow(clippy::useless_conversion)]

use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyOverflowError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString, PyTuple};

use crate::dna::DNA;
use crate::protein::Protein;
use crate::seq_shared;
use crate::utils::{self, PyRnaNeedle};
use biorust_core::seq::rna::RnaSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
pub struct RNA {
    pub(crate) inner: RnaSeq,
}

#[pymethods]
impl RNA {
    #[new]
    fn new(seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let bytes: Vec<u8> = if let Ok(s) = seq.downcast::<PyString>() {
            s.to_str()?.as_bytes().to_vec()
        } else {
            seq.extract::<Vec<u8>>()
                .map_err(|_| PyValueError::new_err("RNA() expects str or bytes-like input"))?
        };

        let inner = RnaSeq::new(bytes).map_err(|e| PyValueError::new_err(e.to_string()))?;
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

    fn back_transcribe(&self) -> DNA {
        DNA {
            inner: self.inner.back_transcribe(),
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
        seq_shared::seq_to_bytes(py, self.as_bytes())
    }

    fn __len__(&self) -> usize {
        self.as_bytes().len()
    }

    fn __richcmp__(&self, other: PyRef<'_, RNA>, op: CompareOp) -> PyResult<bool> {
        let other = other.as_bytes();
        match op {
            CompareOp::Eq => Ok(self.as_bytes() == other),
            CompareOp::Ne => Ok(self.as_bytes() != other),
            CompareOp::Lt => Ok(self.as_bytes() < other),
            CompareOp::Le => Ok(self.as_bytes() <= other),
            CompareOp::Gt => Ok(self.as_bytes() > other),
            CompareOp::Ge => Ok(self.as_bytes() >= other),
        }
    }

    fn __bytes__<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        seq_shared::seq_to_bytes(py, self.as_bytes())
    }

    fn __str__(&self) -> PyResult<String> {
        seq_shared::seq_str(self.as_bytes())
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(seq_shared::seq_repr(self.as_bytes(), "RNA"))
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> RNAIterator {
        RNAIterator {
            bytes: slf.as_bytes().to_vec(),
            index: 0,
        }
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let make = |out: Vec<u8>| -> PyResult<PyObject> {
            let inner = RnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Py::new(py, RNA { inner })?.to_object(py))
        };

        seq_shared::seq_getitem(self.as_bytes(), index, make)
    }

    fn __add__(&self, other: PyRef<'_, RNA>) -> PyResult<Self> {
        let inner = concat_rna_bytes(self.as_bytes(), other.as_bytes())?;
        Ok(Self { inner })
    }

    fn __mul__(&self, num: isize) -> PyResult<Self> {
        let inner = repeat_rna_bytes(self.as_bytes(), num)?;
        Ok(Self { inner })
    }

    fn __rmul__(&self, num: isize) -> PyResult<Self> {
        self.__mul__(num)
    }

    fn count(&self, sub: &Bound<'_, PyAny>) -> PyResult<usize> {
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.count(&other.inner),
            PyRnaNeedle::Bytes(bytes) => self.inner.count(bytes.as_slice()),
            PyRnaNeedle::Byte(b) => self.inner.count(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn count_overlap(&self, sub: &Bound<'_, PyAny>) -> PyResult<usize> {
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.count_overlap(&other.inner),
            PyRnaNeedle::Bytes(bytes) => self.inner.count_overlap(bytes.as_slice()),
            PyRnaNeedle::Byte(b) => self.inner.count_overlap(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __contains__(&self, sub: &Bound<'_, PyAny>) -> PyResult<bool> {
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.contains(&other.inner),
            PyRnaNeedle::Bytes(bytes) => self.inner.contains(bytes.as_slice()),
            PyRnaNeedle::Byte(b) => self.inner.contains(b),
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
        let window = seq_shared::startswith_window(self.as_bytes(), start, end);
        let matches = |needle: PyRnaNeedle<'_>| -> bool {
            let needle = rna_needle_bytes(&needle);
            seq_shared::needle_starts_with(window, needle)
        };

        if let Ok(tuple) = prefix.downcast::<PyTuple>() {
            for item in tuple.iter() {
                let needle = utils::extract_rna_needle(&item)?;
                if matches(needle) {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        let needle = utils::extract_rna_needle(prefix)?;
        Ok(matches(needle))
    }

    #[pyo3(signature = (suffix, start=None, end=None))]
    fn endswith(
        &self,
        suffix: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<bool> {
        let window = seq_shared::startswith_window(self.as_bytes(), start, end);
        let matches = |needle: PyRnaNeedle<'_>| -> bool {
            let needle = rna_needle_bytes(&needle);
            seq_shared::needle_ends_with(window, needle)
        };

        if let Ok(tuple) = suffix.downcast::<PyTuple>() {
            for item in tuple.iter() {
                let needle = utils::extract_rna_needle(&item)?;
                if matches(needle) {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        let needle = utils::extract_rna_needle(suffix)?;
        Ok(matches(needle))
    }

    #[pyo3(signature = (sep=None, maxsplit=-1))]
    fn split<'py>(
        &self,
        py: Python<'py>,
        sep: Option<&Bound<'py, PyAny>>,
        maxsplit: isize,
    ) -> PyResult<Vec<Py<RNA>>> {
        let bytes = self.as_bytes();
        let parts = match sep {
            None => seq_shared::split_on_whitespace(bytes, maxsplit),
            Some(obj) => {
                let needle = utils::extract_rna_needle(obj)?;
                seq_shared::split_on_sep(bytes, rna_needle_bytes(&needle), maxsplit)?
            }
        };

        seq_shared::list_from_parts(parts, |part| {
            let inner = RnaSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Py::new(py, RNA { inner })
        })
    }

    #[pyo3(signature = (sep=None, maxsplit=-1))]
    fn rsplit<'py>(
        &self,
        py: Python<'py>,
        sep: Option<&Bound<'py, PyAny>>,
        maxsplit: isize,
    ) -> PyResult<Vec<Py<RNA>>> {
        let bytes = self.as_bytes();
        let parts = match sep {
            None => seq_shared::rsplit_on_whitespace(bytes, maxsplit),
            Some(obj) => {
                let needle = utils::extract_rna_needle(obj)?;
                seq_shared::rsplit_on_sep(bytes, rna_needle_bytes(&needle), maxsplit)?
            }
        };

        seq_shared::list_from_parts(parts, |part| {
            let inner = RnaSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Py::new(py, RNA { inner })
        })
    }

    #[pyo3(signature = (chars=None))]
    fn strip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_rna_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(rna_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, true, true)?;
        let inner = RnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn lstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_rna_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(rna_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, true, false)?;
        let inner = RnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn rstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_rna_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(rna_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, false, true)?;
        let inner = RnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn upper(&self) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = RnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };
        seq_shared::seq_upper(self.as_bytes(), make)
    }

    fn lower(&self) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = RnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };
        seq_shared::seq_lower(self.as_bytes(), make)
    }

    #[pyo3(signature = (sub, start=None, end=None))]
    fn find(
        &self,
        sub: &Bound<'_, PyAny>,
        start: Option<isize>,
        end: Option<isize>,
    ) -> PyResult<isize> {
        let (s, e) = utils::normalize_range(self.as_bytes().len(), start, end);
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.find(&other.inner, s, e),
            PyRnaNeedle::Bytes(bytes) => self.inner.find(bytes.as_slice(), s, e),
            PyRnaNeedle::Byte(b) => self.inner.find(b, s, e),
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
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.find(&other.inner, s, e),
            PyRnaNeedle::Bytes(bytes) => self.inner.find(bytes.as_slice(), s, e),
            PyRnaNeedle::Byte(b) => self.inner.find(b, s, e),
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
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.rfind(&other.inner, s, e),
            PyRnaNeedle::Bytes(bytes) => self.inner.rfind(bytes.as_slice(), s, e),
            PyRnaNeedle::Byte(b) => self.inner.rfind(b, s, e),
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
        let needle = utils::extract_rna_needle(sub)?;

        let res = match needle {
            PyRnaNeedle::Rna(other) => self.inner.rfind(&other.inner, s, e),
            PyRnaNeedle::Bytes(bytes) => self.inner.rfind(bytes.as_slice(), s, e),
            PyRnaNeedle::Byte(b) => self.inner.rfind(b, s, e),
        };

        match res.map_err(|err| PyValueError::new_err(err.to_string()))? {
            Some(pos) => Ok(pos as isize),
            None => Err(PyValueError::new_err("subsection not found")),
        }
    }
}

#[pyclass]
struct RNAIterator {
    bytes: Vec<u8>,
    index: usize,
}

#[pymethods]
impl RNAIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<String> {
        if self.index < self.bytes.len() {
            let ch = self.bytes[self.index] as char;
            self.index += 1;
            Some(ch.to_string())
        } else {
            None
        }
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RNA>()?;
    m.add_class::<RNAIterator>()?;
    Ok(())
}

fn rna_needle_bytes<'a>(needle: &'a PyRnaNeedle<'a>) -> seq_shared::NeedleBytes<'a> {
    match needle {
        PyRnaNeedle::Rna(other) => seq_shared::NeedleBytes::Bytes(other.as_bytes()),
        PyRnaNeedle::Bytes(bytes) => seq_shared::NeedleBytes::Bytes(bytes.as_slice()),
        PyRnaNeedle::Byte(b) => seq_shared::NeedleBytes::Byte(*b),
    }
}

fn concat_rna_bytes(left: &[u8], right: &[u8]) -> PyResult<RnaSeq> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    out.extend_from_slice(left);
    out.extend_from_slice(right);
    RnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))
}

fn repeat_rna_bytes(bytes: &[u8], n: isize) -> PyResult<RnaSeq> {
    if n <= 0 || bytes.is_empty() {
        return RnaSeq::new(Vec::new()).map_err(|e| PyValueError::new_err(e.to_string()));
    }

    let n = n as usize;
    let total = bytes
        .len()
        .checked_mul(n)
        .ok_or_else(|| PyOverflowError::new_err("repeat count causes overflow"))?;
    let mut out = Vec::with_capacity(total);
    for _ in 0..n {
        out.extend_from_slice(bytes);
    }
    RnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))
}
