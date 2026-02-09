#![allow(clippy::useless_conversion)]

use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyOverflowError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString, PyTuple};

use crate::seq_shared;
use crate::utils::{self, PyProteinNeedle};
use biorust_core::seq::protein::ProteinSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
#[derive(Clone)]
pub struct Protein {
    pub(crate) inner: ProteinSeq,
}

#[pymethods]
impl Protein {
    #[new]
    fn new(seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let bytes: Vec<u8> = if let Ok(s) = seq.downcast::<PyString>() {
            s.to_str()?.as_bytes().to_vec()
        } else {
            seq.extract::<Vec<u8>>()
                .map_err(|_| PyValueError::new_err("Protein() expects str or bytes-like input"))?
        };

        let inner = ProteinSeq::new(bytes).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
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

    fn __richcmp__(&self, other: PyRef<'_, Protein>, op: CompareOp) -> PyResult<bool> {
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
        Ok(seq_shared::seq_repr(self.as_bytes(), "Protein"))
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<ProteinIterator> {
        let py = slf.py();
        Ok(ProteinIterator {
            seq: Py::new(py, slf.clone())?,
            index: 0,
        })
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let make = |out: Vec<u8>| -> PyResult<PyObject> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Py::new(py, Protein { inner })?.to_object(py))
        };

        seq_shared::seq_getitem(self.as_bytes(), index, make)
    }

    fn __add__(&self, other: PyRef<'_, Protein>) -> PyResult<Self> {
        let inner = concat_protein_bytes(self.as_bytes(), other.as_bytes())?;
        Ok(Self { inner })
    }

    fn __mul__(&self, num: isize) -> PyResult<Self> {
        let inner = repeat_protein_bytes(self.as_bytes(), num)?;
        Ok(Self { inner })
    }

    fn __rmul__(&self, num: isize) -> PyResult<Self> {
        self.__mul__(num)
    }

    fn count(&self, sub: &Bound<'_, PyAny>) -> PyResult<usize> {
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.count(&other.inner),
            PyProteinNeedle::Bytes(bytes) => self.inner.count(bytes.as_slice()),
            PyProteinNeedle::Byte(b) => self.inner.count(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn count_overlap(&self, sub: &Bound<'_, PyAny>) -> PyResult<usize> {
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.count_overlap(&other.inner),
            PyProteinNeedle::Bytes(bytes) => self.inner.count_overlap(bytes.as_slice()),
            PyProteinNeedle::Byte(b) => self.inner.count_overlap(b),
        };

        res.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __contains__(&self, sub: &Bound<'_, PyAny>) -> PyResult<bool> {
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.contains(&other.inner),
            PyProteinNeedle::Bytes(bytes) => self.inner.contains(bytes.as_slice()),
            PyProteinNeedle::Byte(b) => self.inner.contains(b),
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
        let matches = |needle: PyProteinNeedle<'_>| -> bool {
            let needle = protein_needle_bytes(&needle);
            seq_shared::needle_starts_with(window, needle)
        };

        if let Ok(tuple) = prefix.downcast::<PyTuple>() {
            for item in tuple.iter() {
                let needle = utils::extract_protein_needle(&item)?;
                if matches(needle) {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        let needle = utils::extract_protein_needle(prefix)?;
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
        let matches = |needle: PyProteinNeedle<'_>| -> bool {
            let needle = protein_needle_bytes(&needle);
            seq_shared::needle_ends_with(window, needle)
        };

        if let Ok(tuple) = suffix.downcast::<PyTuple>() {
            for item in tuple.iter() {
                let needle = utils::extract_protein_needle(&item)?;
                if matches(needle) {
                    return Ok(true);
                }
            }
            return Ok(false);
        }

        let needle = utils::extract_protein_needle(suffix)?;
        Ok(matches(needle))
    }

    #[pyo3(signature = (sep=None, maxsplit=-1))]
    fn split<'py>(
        &self,
        py: Python<'py>,
        sep: Option<&Bound<'py, PyAny>>,
        maxsplit: isize,
    ) -> PyResult<Vec<Py<Protein>>> {
        let bytes = self.as_bytes();
        let parts = match sep {
            None => seq_shared::split_on_whitespace(bytes, maxsplit),
            Some(obj) => {
                let needle = utils::extract_protein_needle(obj)?;
                seq_shared::split_on_sep(bytes, protein_needle_bytes(&needle), maxsplit)?
            }
        };

        seq_shared::list_from_parts(parts, |part| {
            let inner = ProteinSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Py::new(py, Protein { inner })
        })
    }

    #[pyo3(signature = (sep=None, maxsplit=-1))]
    fn rsplit<'py>(
        &self,
        py: Python<'py>,
        sep: Option<&Bound<'py, PyAny>>,
        maxsplit: isize,
    ) -> PyResult<Vec<Py<Protein>>> {
        let bytes = self.as_bytes();
        let parts = match sep {
            None => seq_shared::rsplit_on_whitespace(bytes, maxsplit),
            Some(obj) => {
                let needle = utils::extract_protein_needle(obj)?;
                seq_shared::rsplit_on_sep(bytes, protein_needle_bytes(&needle), maxsplit)?
            }
        };

        seq_shared::list_from_parts(parts, |part| {
            let inner = ProteinSeq::new(part).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Py::new(py, Protein { inner })
        })
    }

    #[pyo3(signature = (chars=None))]
    fn strip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_protein_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(protein_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, true, true)?;
        let inner = ProteinSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn lstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_protein_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(protein_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, true, false)?;
        let inner = ProteinSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[pyo3(signature = (chars=None))]
    fn rstrip(&self, chars: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let bytes = self.as_bytes();
        let needle = match chars {
            Some(obj) => Some(utils::extract_protein_needle(obj)?),
            None => None,
        };
        let needle = needle.as_ref().map(protein_needle_bytes);
        let (start, end) = seq_shared::trim_range(bytes, needle, false, true)?;
        let inner = ProteinSeq::new(bytes[start..end].to_vec())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn upper(&self) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };
        seq_shared::seq_upper(self.as_bytes(), make)
    }

    fn lower(&self) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
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
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.find(&other.inner, s, e),
            PyProteinNeedle::Bytes(bytes) => self.inner.find(bytes.as_slice(), s, e),
            PyProteinNeedle::Byte(b) => self.inner.find(b, s, e),
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
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.find(&other.inner, s, e),
            PyProteinNeedle::Bytes(bytes) => self.inner.find(bytes.as_slice(), s, e),
            PyProteinNeedle::Byte(b) => self.inner.find(b, s, e),
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
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.rfind(&other.inner, s, e),
            PyProteinNeedle::Bytes(bytes) => self.inner.rfind(bytes.as_slice(), s, e),
            PyProteinNeedle::Byte(b) => self.inner.rfind(b, s, e),
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
        let needle = utils::extract_protein_needle(sub)?;

        let res = match needle {
            PyProteinNeedle::Protein(other) => self.inner.rfind(&other.inner, s, e),
            PyProteinNeedle::Bytes(bytes) => self.inner.rfind(bytes.as_slice(), s, e),
            PyProteinNeedle::Byte(b) => self.inner.rfind(b, s, e),
        };

        match res.map_err(|err| PyValueError::new_err(err.to_string()))? {
            Some(pos) => Ok(pos as isize),
            None => Err(PyValueError::new_err("subsection not found")),
        }
    }

    fn reverse(&self) -> Self {
        Self {
            inner: self.inner.reverse(),
        }
    }

    fn counts(&self) -> Vec<(String, u32)> {
        let counts = self.inner.counts();
        counts
            .iter()
            .enumerate()
            .filter_map(|(b, &count)| {
                if count == 0 {
                    None
                } else {
                    Some(((b as u8 as char).to_string(), count))
                }
            })
            .collect()
    }

    fn frequencies(&self) -> Vec<(String, f64)> {
        let freq = self.inner.frequencies();
        freq.iter()
            .enumerate()
            .filter_map(|(b, &val)| {
                if val == 0.0 {
                    None
                } else {
                    Some(((b as u8 as char).to_string(), val))
                }
            })
            .collect()
    }

    fn aa_counts_20(&self) -> Vec<(String, u32)> {
        let counts = self.inner.aa_counts_20();
        let letters = b"ARNDCEQGHILKMFPSTWYV";
        letters
            .iter()
            .enumerate()
            .map(|(i, &b)| ((b as char).to_string(), counts[i]))
            .collect()
    }

    fn aa_frequencies_20(&self) -> Vec<(String, f64)> {
        let freq = self.inner.aa_frequencies_20();
        let letters = b"ARNDCEQGHILKMFPSTWYV";
        letters
            .iter()
            .enumerate()
            .map(|(i, &b)| ((b as char).to_string(), freq[i]))
            .collect()
    }

    fn shannon_entropy(&self) -> f64 {
        self.inner.shannon_entropy()
    }

    fn molecular_weight(&self) -> PyResult<f64> {
        self.inner
            .molecular_weight()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn hydrophobicity(&self) -> PyResult<f64> {
        self.inner
            .hydrophobicity()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn hydrophobicity_profile(&self, window: usize) -> PyResult<Vec<f64>> {
        self.inner
            .hydrophobicity_profile(window)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn net_charge(&self, ph: f64) -> PyResult<f64> {
        self.inner
            .net_charge(ph)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn isoelectric_point(&self) -> PyResult<f64> {
        self.inner
            .isoelectric_point()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn validate_strict_20(&self) -> PyResult<()> {
        self.inner
            .validate_strict_20()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn has_ambiguous(&self) -> bool {
        self.inner.has_ambiguous()
    }

    fn unknown_positions(&self) -> Vec<usize> {
        self.inner.unknown_positions()
    }
}

#[pyclass]
struct ProteinIterator {
    seq: Py<Protein>,
    index: usize,
}

#[pymethods]
impl ProteinIterator {
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

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Protein>()?;
    m.add_class::<ProteinIterator>()?;
    Ok(())
}

fn protein_needle_bytes<'a>(needle: &'a PyProteinNeedle<'a>) -> seq_shared::NeedleBytes<'a> {
    match needle {
        PyProteinNeedle::Protein(other) => seq_shared::NeedleBytes::Bytes(other.as_bytes()),
        PyProteinNeedle::Bytes(bytes) => seq_shared::NeedleBytes::Bytes(bytes.as_slice()),
        PyProteinNeedle::Byte(b) => seq_shared::NeedleBytes::Byte(*b),
    }
}

fn concat_protein_bytes(left: &[u8], right: &[u8]) -> PyResult<ProteinSeq> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    out.extend_from_slice(left);
    out.extend_from_slice(right);
    ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))
}

fn repeat_protein_bytes(bytes: &[u8], n: isize) -> PyResult<ProteinSeq> {
    if n <= 0 || bytes.is_empty() {
        return ProteinSeq::new(Vec::new()).map_err(|e| PyValueError::new_err(e.to_string()));
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
    ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))
}
