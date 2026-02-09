#![allow(clippy::useless_conversion)]

use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyOverflowError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString, PyTuple};

use crate::protein::Protein;
use crate::rna::RNA;
use crate::seq_shared;
use crate::utils::{self, PyDnaNeedle};
use biorust_core::seq::dna::DnaSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
#[derive(Clone)]
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

    fn transcribe(&self) -> RNA {
        RNA {
            inner: self.inner.transcribe(),
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

    fn __richcmp__(&self, other: PyRef<'_, DNA>, op: CompareOp) -> PyResult<bool> {
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
        Ok(seq_shared::seq_repr(self.as_bytes(), "DNA"))
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<DNAIterator> {
        let py = slf.py();
        Ok(DNAIterator {
            seq: Py::new(py, slf.clone())?,
            index: 0,
        })
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let make = |out: Vec<u8>| -> PyResult<PyObject> {
            let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Py::new(py, DNA { inner })?.to_object(py))
        };

        seq_shared::seq_getitem(self.as_bytes(), index, make)
    }

    fn __add__(&self, other: PyRef<'_, DNA>) -> PyResult<Self> {
        let inner = concat_dna_bytes(self.as_bytes(), other.as_bytes())?;
        Ok(Self { inner })
    }

    fn __mul__(&self, num: isize) -> PyResult<Self> {
        let inner = repeat_dna_bytes(self.as_bytes(), num)?;
        Ok(Self { inner })
    }

    fn __rmul__(&self, num: isize) -> PyResult<Self> {
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
        let window = seq_shared::startswith_window(self.as_bytes(), start, end);
        let matches = |needle: PyDnaNeedle<'_>| -> bool {
            let needle = dna_needle_bytes(&needle);
            seq_shared::needle_starts_with(window, needle)
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
        let window = seq_shared::startswith_window(self.as_bytes(), start, end);
        let matches = |needle: PyDnaNeedle<'_>| -> bool {
            let needle = dna_needle_bytes(&needle);
            seq_shared::needle_ends_with(window, needle)
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
            None => seq_shared::split_on_whitespace(bytes, maxsplit),
            Some(obj) => {
                let needle = utils::extract_dna_needle(obj)?;
                seq_shared::split_on_sep(bytes, dna_needle_bytes(&needle), maxsplit)?
            }
        };

        seq_shared::list_from_parts(parts, |part| {
            let inner = DnaSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Py::new(py, DNA { inner })
        })
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
            None => seq_shared::rsplit_on_whitespace(bytes, maxsplit),
            Some(obj) => {
                let needle = utils::extract_dna_needle(obj)?;
                seq_shared::rsplit_on_sep(bytes, dna_needle_bytes(&needle), maxsplit)?
            }
        };

        seq_shared::list_from_parts(parts, |part| {
            let inner = DnaSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Py::new(py, DNA { inner })
        })
    }

    #[pyo3(signature = (chars=None))]
    fn strip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_dna_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(dna_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, true, true)?;
        let inner = DnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn lstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_dna_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(dna_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, true, false)?;
        let inner = DnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn rstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_dna_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(dna_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, false, true)?;
        let inner = DnaSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn upper(&self) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };
        seq_shared::seq_upper(self.as_bytes(), make)
    }

    fn lower(&self) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
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

#[pyclass]
struct DNAIterator {
    seq: Py<DNA>,
    index: usize,
}

#[pymethods]
impl DNAIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<String> {
        let borrowed = self.seq.borrow(py);
        let bytes = borrowed.as_bytes();
        if self.index < bytes.len() {
            let ch = bytes[self.index] as char;
            self.index += 1;
            Some(ch.to_string())
        } else {
            None
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

    let bytes = utils::extract_dna_bytes(seq)?;
    let inner = DnaSeq::new(bytes).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(DNA {
        inner: inner.complement(),
    })
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DNA>()?;
    m.add_class::<DNAIterator>()?;
    m.add_function(wrap_pyfunction!(complement, m)?)?;
    Ok(())
}

fn dna_needle_bytes<'a>(needle: &'a PyDnaNeedle<'a>) -> seq_shared::NeedleBytes<'a> {
    match needle {
        PyDnaNeedle::Dna(other) => seq_shared::NeedleBytes::Bytes(other.as_bytes()),
        PyDnaNeedle::Bytes(bytes) => seq_shared::NeedleBytes::Bytes(bytes.as_slice()),
        PyDnaNeedle::Byte(b) => seq_shared::NeedleBytes::Byte(*b),
    }
}

fn concat_dna_bytes(left: &[u8], right: &[u8]) -> PyResult<DnaSeq> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    out.extend_from_slice(left);
    out.extend_from_slice(right);
    DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))
}

fn repeat_dna_bytes(bytes: &[u8], n: isize) -> PyResult<DnaSeq> {
    if n <= 0 || bytes.is_empty() {
        return DnaSeq::new(Vec::new()).map_err(|e| PyValueError::new_err(e.to_string()));
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
    DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))
}
