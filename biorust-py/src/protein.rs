#![allow(clippy::useless_conversion)]

use pyo3::basic::CompareOp;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString, PyTuple};

use crate::seq_shared;
use crate::utils::{self, PyProteinNeedle};
use biorust_core::seq::protein::ProteinSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
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

    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: CompareOp) -> PyResult<bool> {
        let other = utils::extract_protein_bytes(other)?;

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
        seq_shared::seq_to_bytes(py, self.as_bytes())
    }

    fn __str__(&self) -> PyResult<String> {
        seq_shared::seq_str(self.as_bytes())
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(seq_shared::seq_repr(self.as_bytes(), "Protein"))
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let make = |out: Vec<u8>| -> PyResult<PyObject> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Py::new(py, Protein { inner })?.to_object(py))
        };

        seq_shared::seq_getitem(self.as_bytes(), index, make)
    }

    fn __radd__(&self, seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };

        seq_shared::seq_radd(seq, self.as_bytes(), utils::extract_protein_bytes, make)
    }

    fn __add__(&self, seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };

        seq_shared::seq_add(self.as_bytes(), seq, utils::extract_protein_bytes, make)
    }

    fn __mul__(&self, num: &Bound<'_, PyAny>) -> PyResult<Self> {
        let make = |out: Vec<u8>| -> PyResult<Self> {
            let inner = ProteinSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
            Ok(Self { inner })
        };

        seq_shared::seq_mul(self.as_bytes(), num, make)
    }

    fn __rmul__(&self, num: &Bound<'_, PyAny>) -> PyResult<Self> {
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
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Protein>()?;
    Ok(())
}

fn protein_needle_bytes<'a>(needle: &'a PyProteinNeedle<'a>) -> seq_shared::NeedleBytes<'a> {
    match needle {
        PyProteinNeedle::Protein(other) => seq_shared::NeedleBytes::Bytes(other.as_bytes()),
        PyProteinNeedle::Bytes(bytes) => seq_shared::NeedleBytes::Bytes(bytes.as_slice()),
        PyProteinNeedle::Byte(b) => seq_shared::NeedleBytes::Byte(*b),
    }
}
