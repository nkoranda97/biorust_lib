use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString};

use crate::protein::Protein;
use crate::seq_shared;
use biorust_core::seq::gapped_protein::GappedProteinSeq;

#[pyclass(frozen)]
#[derive(Clone)]
pub struct GappedProtein {
    pub(crate) inner: GappedProteinSeq,
}

#[pymethods]
impl GappedProtein {
    #[new]
    fn new(seq: &Bound<'_, PyAny>) -> PyResult<Self> {
        let bytes: Vec<u8> = if let Ok(s) = seq.downcast::<PyString>() {
            s.to_str()?.as_bytes().to_vec()
        } else {
            seq.extract::<Vec<u8>>().map_err(|_| {
                pyo3::exceptions::PyValueError::new_err(
                    "GappedProtein() expects str or bytes-like input",
                )
            })?
        };

        let inner = GappedProteinSeq::new(bytes)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    fn ungapped(&self) -> Protein {
        Protein {
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

    fn __richcmp__(&self, other: PyRef<'_, GappedProtein>, op: CompareOp) -> PyResult<bool> {
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
        Ok(seq_shared::seq_repr(self.as_bytes(), "GappedProtein"))
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let make = |out: Vec<u8>| -> PyResult<PyObject> {
            let inner = GappedProteinSeq::new(out)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
            Ok(Py::new(py, GappedProtein { inner })?.to_object(py))
        };

        seq_shared::seq_getitem(self.as_bytes(), index, make)
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GappedProtein>()?;
    Ok(())
}
