use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString};

use crate::dna::DNA;
use crate::seq_shared;
use biorust_core::seq::gapped_dna::GappedDnaSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
#[derive(Clone)]
pub struct GappedDNA {
    pub(crate) inner: GappedDnaSeq,
}

#[pymethods]
impl GappedDNA {
    #[new]
    fn new(seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let bytes: Vec<u8> = if let Ok(s) = seq.downcast::<PyString>() {
            s.to_str()?.as_bytes().to_vec()
        } else {
            seq.extract::<Vec<u8>>().map_err(|_| {
                pyo3::exceptions::PyValueError::new_err(
                    "GappedDNA() expects str or bytes-like input",
                )
            })?
        };

        let inner = GappedDnaSeq::new(bytes)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn ungapped(&self) -> DNA {
        DNA {
            inner: self.inner.ungapped(),
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

    fn __richcmp__(&self, other: PyRef<'_, GappedDNA>, op: CompareOp) -> PyResult<bool> {
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
        Ok(seq_shared::seq_repr(self.as_bytes(), "GappedDNA"))
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let make = |out: Vec<u8>| -> PyResult<PyObject> {
            let inner = GappedDnaSeq::new(out)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
            Ok(Py::new(py, GappedDNA { inner })?.to_object(py))
        };

        seq_shared::seq_getitem(self.as_bytes(), index, make)
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GappedDNA>()?;
    Ok(())
}
