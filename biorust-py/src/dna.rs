#![allow(clippy::useless_conversion)]

use pyo3::basic::CompareOp;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PySlice, PyString};

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
}

#[pyfunction]
fn complement(a: u8) -> u8 {
    biorust_core::alphabets::dna::complement(a)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DNA>()?;
    m.add_function(wrap_pyfunction!(complement, m)?)?;
    Ok(())
}
