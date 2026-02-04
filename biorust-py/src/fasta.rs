#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::dna_record_batch::DNARecordBatch;
use biorust_core::error::BioError;
use biorust_core::io::fasta;
use biorust_core::seq::dna::DnaSeq;

#[pyfunction]
fn read_fasta(py: Python<'_>, path: &str) -> PyResult<PyObject> {
    let batch = fasta::read_fasta_batch_from_path::<DnaSeq>(path).map_err(|err| match err {
        BioError::FastaIo(io) => PyIOError::new_err(io.to_string()),
        other => PyValueError::new_err(other.to_string()),
    })?;
    let out = DNARecordBatch { inner: batch };
    Ok(Py::new(py, out)?.to_object(py))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    Ok(())
}
