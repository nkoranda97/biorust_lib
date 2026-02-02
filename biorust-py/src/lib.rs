use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

use biorust_core::seq::dna::DnaSeq;

fn extract_bytes<'py>(seq: &Bound<'py, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(dna) = seq.extract::<PyRef<'py, DNA>>() {
        return Ok(dna.inner.as_bytes().to_vec());
    }
    if let Ok(s) = seq.downcast::<PyString>() {
        return Ok(s.to_str()?.as_bytes().to_vec());
    }
    seq.extract::<Vec<u8>>()
        .map_err(|_| PyTypeError::new_err("expected DNA, str, or bytes-like object"))
}

#[pyclass(frozen)]
pub struct DNA {
    inner: DnaSeq, // <-- store core struct, not Vec<u8>
}

#[pymethods]
impl DNA {
    #[new]
    fn new<'py>(seq: &Bound<'py, PyAny>) -> PyResult<Self> {
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

    fn to_bytes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, self.inner.as_bytes())
    }

    fn __len__(&self) -> usize {
        self.inner.as_bytes().len()
    }

    fn __bytes__<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new_bound(py, self.inner.as_bytes())
    }

    fn __str__(&self) -> PyResult<String> {
        std::str::from_utf8(self.inner.as_bytes())
            .map(|s| s.to_string())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> PyResult<String> {
        let s = std::str::from_utf8(self.inner.as_bytes()).unwrap_or("<bytes>");
        Ok(format!("DNA({s:?})"))
    }

    fn __radd__<'py>(&self, seq: &Bound<'py, PyAny>) -> PyResult<Self> {
        let mut left = extract_bytes(seq)?;

        // left + self
        left.extend_from_slice(self.inner.as_bytes());

        let inner = DnaSeq::new(left).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn __add__<'py>(&self, seq: &Bound<'py, PyAny>) -> PyResult<Self> {
        let right = extract_bytes(seq)?;

        // self + right
        let mut out = Vec::with_capacity(self.inner.as_bytes().len() + right.len());
        out.extend_from_slice(self.inner.as_bytes());
        out.extend_from_slice(&right);

        let inner = DnaSeq::new(out).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }
}

#[pyfunction]
fn complement(a: u8) -> u8 {
    biorust_core::alphabets::dna::complement(a)
}

#[pymodule]
fn _native(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DNA>()?;
    m.add_function(wrap_pyfunction!(complement, m)?)?;
    Ok(())
}