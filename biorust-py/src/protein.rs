use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyModule, PyString};

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
        PyBytes::new_bound(py, self.as_bytes())
    }

    fn __len__(&self) -> usize {
        self.as_bytes().len()
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
        Ok(format!("Protein({s:?})"))
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Protein>()?;
    Ok(())
}
