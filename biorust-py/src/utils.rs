use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::dna::DNA;

pub fn extract_bytes<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(dna) = obj.extract::<PyRef<'py, DNA>>() {
        return Ok(dna.inner.as_bytes().to_vec());
    }

    if let Ok(s) = obj.downcast::<PyString>() {
        return Ok(s.to_str()?.as_bytes().to_vec());
    }

    obj.extract::<Vec<u8>>()
        .map_err(|_| PyTypeError::new_err("expected DNA, str, or bytes-like object"))
}
