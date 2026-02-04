use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::dna::DNA;
use crate::protein::Protein;

pub fn extract_dna_bytes<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Vec<u8>> {
    if let Ok(dna) = obj.extract::<PyRef<'py, DNA>>() {
        return Ok(dna.as_bytes().to_vec());
    }

    if let Ok(s) = obj.downcast::<PyString>() {
        return Ok(s.to_str()?.as_bytes().to_vec());
    }

    obj.extract::<Vec<u8>>()
        .map_err(|_| PyTypeError::new_err("expected DNA, str, or bytes-like object"))
}

pub enum PyDnaNeedle<'a> {
    Bytes(Vec<u8>),
    Byte(u8),
    Dna(PyRef<'a, DNA>),
}

pub enum PyProteinNeedle<'a> {
    Bytes(Vec<u8>),
    Byte(u8),
    Protein(PyRef<'a, Protein>),
}

pub fn extract_dna_needle<'py>(obj: &'py Bound<'py, PyAny>) -> PyResult<PyDnaNeedle<'py>> {
    if let Ok(dna) = obj.extract::<PyRef<'py, DNA>>() {
        return Ok(PyDnaNeedle::Dna(dna));
    }

    if let Ok(s) = obj.downcast::<PyString>() {
        return Ok(PyDnaNeedle::Bytes(s.to_str()?.as_bytes().to_vec()));
    }

    if let Ok(b) = obj.extract::<Vec<u8>>() {
        return Ok(PyDnaNeedle::Bytes(b));
    }

    if let Ok(n) = obj.extract::<i64>() {
        if (0..=255).contains(&n) {
            return Ok(PyDnaNeedle::Byte(n as u8));
        }
        return Err(PyValueError::new_err("sub must be an int in range 0..=255"));
    }

    Err(PyTypeError::new_err(
        "sub must be DNA, str, bytes-like, or int (0..=255)",
    ))
}

pub fn extract_protein_needle<'py>(obj: &'py Bound<'py, PyAny>) -> PyResult<PyProteinNeedle<'py>> {
    if let Ok(protein) = obj.extract::<PyRef<'py, Protein>>() {
        return Ok(PyProteinNeedle::Protein(protein));
    }

    if let Ok(s) = obj.downcast::<PyString>() {
        return Ok(PyProteinNeedle::Bytes(s.to_str()?.as_bytes().to_vec()));
    }

    if let Ok(b) = obj.extract::<Vec<u8>>() {
        return Ok(PyProteinNeedle::Bytes(b));
    }

    if let Ok(n) = obj.extract::<i64>() {
        if (0..=255).contains(&n) {
            return Ok(PyProteinNeedle::Byte(n as u8));
        }
        return Err(PyValueError::new_err("sub must be an int in range 0..=255"));
    }

    Err(PyTypeError::new_err(
        "sub must be Protein, str, bytes-like, or int (0..=255)",
    ))
}

pub fn normalize_range(len: usize, start: Option<isize>, end: Option<isize>) -> (usize, usize) {
    let n = len as isize;

    let mut s = start.unwrap_or(0);
    let mut e = end.unwrap_or(n);

    if s < 0 {
        s += n;
    }
    if e < 0 {
        e += n;
    }

    s = s.clamp(0, n);
    e = e.clamp(0, n);

    (s as usize, e as usize)
}
