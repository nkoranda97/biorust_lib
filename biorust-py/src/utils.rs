use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::dna::DNA;
use crate::protein::Protein;
use crate::rna::RNA;

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

pub enum PyRnaNeedle<'a> {
    Bytes(Vec<u8>),
    Byte(u8),
    Rna(PyRef<'a, RNA>),
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

pub fn extract_rna_needle<'py>(obj: &'py Bound<'py, PyAny>) -> PyResult<PyRnaNeedle<'py>> {
    if let Ok(rna) = obj.extract::<PyRef<'py, RNA>>() {
        return Ok(PyRnaNeedle::Rna(rna));
    }

    if let Ok(s) = obj.downcast::<PyString>() {
        return Ok(PyRnaNeedle::Bytes(s.to_str()?.as_bytes().to_vec()));
    }

    if let Ok(b) = obj.extract::<Vec<u8>>() {
        return Ok(PyRnaNeedle::Bytes(b));
    }

    if let Ok(n) = obj.extract::<i64>() {
        if (0..=255).contains(&n) {
            return Ok(PyRnaNeedle::Byte(n as u8));
        }
        return Err(PyValueError::new_err("sub must be an int in range 0..=255"));
    }

    Err(PyTypeError::new_err(
        "sub must be RNA, str, bytes-like, or int (0..=255)",
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

pub fn parse_frame(obj: &Bound<'_, PyAny>) -> PyResult<biorust_core::seq::TranslationFrame> {
    use biorust_core::seq::TranslationFrame;

    if let Ok(s) = obj.downcast::<PyString>() {
        let s = s.to_str()?;
        if s.eq_ignore_ascii_case("auto") {
            return Ok(TranslationFrame::Auto);
        }
        return Err(PyValueError::new_err("frame must be 1, 2, 3, or \"auto\""));
    }

    if let Ok(n) = obj.extract::<i64>() {
        return match n {
            1 => Ok(TranslationFrame::One),
            2 => Ok(TranslationFrame::Two),
            3 => Ok(TranslationFrame::Three),
            _ => Err(PyValueError::new_err("frame must be 1, 2, 3, or \"auto\"")),
        };
    }

    Err(PyValueError::new_err("frame must be 1, 2, 3, or \"auto\""))
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
