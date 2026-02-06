#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyIOError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule};

use crate::dna_record::DNARecord;
use crate::dna_record_batch::DNARecordBatch;
use crate::protein_record::ProteinRecord;
use crate::protein_record_batch::ProteinRecordBatch;
use biorust_core::error::BioError;
use biorust_core::io::fasta;
use biorust_core::seq::dna::DnaSeq;
use biorust_core::seq::protein::ProteinSeq;
use biorust_core::seq::record::SeqRecord;

#[pyfunction]
fn read_fasta(py: Python<'_>, path: &str) -> PyResult<PyObject> {
    let batch = fasta::read_fasta_batch_from_path::<DnaSeq>(path).map_err(|err| match err {
        BioError::FastaIo(io) => PyIOError::new_err(io.to_string()),
        other => PyValueError::new_err(other.to_string()),
    })?;
    let out = DNARecordBatch {
        inner: batch,
        skipped: Vec::new(),
    };
    Ok(Py::new(py, out)?.to_object(py))
}

#[pyfunction]
#[pyo3(signature = (path, records, *, line_width=60))]
fn write_fasta(path: &str, records: &Bound<'_, PyAny>, line_width: usize) -> PyResult<()> {
    if let Ok(batch) = records.extract::<PyRef<'_, DNARecordBatch>>() {
        return fasta::write_fasta_batch_to_path(path, &batch.inner, line_width)
            .map_err(map_bio_err);
    }
    if let Ok(batch) = records.extract::<PyRef<'_, ProteinRecordBatch>>() {
        return fasta::write_fasta_batch_to_path(path, &batch.inner, line_width)
            .map_err(map_bio_err);
    }

    let mut dna_records: Vec<SeqRecord<DnaSeq>> = Vec::new();
    let mut protein_records: Vec<SeqRecord<ProteinSeq>> = Vec::new();
    for item in records.iter()? {
        let item = item?;
        if let Ok(record) = item.extract::<PyRef<'_, DNARecord>>() {
            if !protein_records.is_empty() {
                return Err(PyTypeError::new_err(
                    "write_fasta records must be all DNARecord or all ProteinRecord",
                ));
            }
            dna_records.push(record.inner.clone());
            continue;
        }
        if let Ok(record) = item.extract::<PyRef<'_, ProteinRecord>>() {
            if !dna_records.is_empty() {
                return Err(PyTypeError::new_err(
                    "write_fasta records must be all DNARecord or all ProteinRecord",
                ));
            }
            protein_records.push(record.inner.clone());
            continue;
        }
        return Err(PyTypeError::new_err(
            "write_fasta expects DNARecordBatch, ProteinRecordBatch, or an iterable of DNARecord or ProteinRecord",
        ));
    }

    if !dna_records.is_empty() {
        return fasta::write_fasta_records_to_path(path, &dna_records, line_width)
            .map_err(map_bio_err);
    }
    if !protein_records.is_empty() {
        return fasta::write_fasta_records_to_path(path, &protein_records, line_width)
            .map_err(map_bio_err);
    }

    Err(PyTypeError::new_err(
        "write_fasta expects a non-empty iterable of DNARecord or ProteinRecord",
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(write_fasta, m)?)?;
    Ok(())
}

fn map_bio_err(err: BioError) -> PyErr {
    match err {
        BioError::FastaIo(io) => PyIOError::new_err(io.to_string()),
        other => PyValueError::new_err(other.to_string()),
    }
}
