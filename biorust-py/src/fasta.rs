#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyIOError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule};

use crate::dna_record::DNARecord;
use crate::dna_record_batch::DNARecordBatch;
use crate::protein_record::ProteinRecord;
use crate::protein_record_batch::ProteinRecordBatch;
use crate::rna_record::RNARecord;
use crate::rna_record_batch::RNARecordBatch;
use biorust_core::error::BioError;
use biorust_core::io::detect::{detect_seq_type, SeqType};
use biorust_core::io::fasta;
use biorust_core::seq::dna::DnaSeq;
use biorust_core::seq::protein::ProteinSeq;
use biorust_core::seq::record::SeqRecord;
use biorust_core::seq::rna::RnaSeq;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[pyfunction]
#[pyo3(signature = (path, *, alphabet="auto"))]
fn read_fasta(py: Python<'_>, path: &str, alphabet: &str) -> PyResult<PyObject> {
    let alpha = match alphabet.to_ascii_lowercase().as_str() {
        "auto" => detect_fasta_type(path)?,
        "dna" => SeqType::Dna,
        "rna" => SeqType::Rna,
        "protein" => SeqType::Protein,
        _ => {
            return Err(PyValueError::new_err(
                "alphabet must be 'auto', 'dna', 'rna', or 'protein'",
            ))
        }
    };

    match alpha {
        SeqType::Dna => {
            let batch = fasta::read_fasta_batch_from_path::<DnaSeq>(path).map_err(map_bio_err)?;
            let out = DNARecordBatch {
                inner: batch,
                skipped: Vec::new(),
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
        SeqType::Rna => {
            let batch = fasta::read_fasta_batch_from_path::<RnaSeq>(path).map_err(map_bio_err)?;
            let out = RNARecordBatch {
                inner: batch,
                skipped: Vec::new(),
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
        SeqType::Protein => {
            let batch =
                fasta::read_fasta_batch_from_path::<ProteinSeq>(path).map_err(map_bio_err)?;
            let out = ProteinRecordBatch {
                inner: batch,
                skipped: Vec::new(),
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
    }
}

/// Peek at the first FASTA record's sequence bytes to detect the alphabet.
fn detect_fasta_type(path: &str) -> PyResult<SeqType> {
    let file = File::open(path).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let reader = BufReader::new(file);
    let mut seq_bytes = Vec::new();
    let mut in_seq = false;

    for line in reader.lines() {
        let line = line.map_err(|e| PyIOError::new_err(e.to_string()))?;
        if line.starts_with('>') {
            if in_seq {
                break; // done with first record
            }
            in_seq = true;
            continue;
        }
        if in_seq {
            for b in line.bytes() {
                if !b.is_ascii_whitespace() {
                    seq_bytes.push(b);
                }
            }
            if seq_bytes.len() >= 1000 {
                break;
            }
        }
    }

    Ok(detect_seq_type(&seq_bytes))
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
    if let Ok(batch) = records.extract::<PyRef<'_, RNARecordBatch>>() {
        return fasta::write_fasta_batch_to_path(path, &batch.inner, line_width)
            .map_err(map_bio_err);
    }

    #[derive(Clone, Copy)]
    enum RecordKind {
        Dna,
        Rna,
        Protein,
    }

    let mut dna_records: Vec<SeqRecord<DnaSeq>> = Vec::new();
    let mut rna_records: Vec<SeqRecord<RnaSeq>> = Vec::new();
    let mut protein_records: Vec<SeqRecord<ProteinSeq>> = Vec::new();
    let mut kind: Option<RecordKind> = None;

    for item in records.iter()? {
        let item = item?;
        if let Ok(record) = item.extract::<PyRef<'_, DNARecord>>() {
            if matches!(kind, Some(RecordKind::Rna | RecordKind::Protein)) {
                return Err(PyTypeError::new_err(
                    "write_fasta records must all be the same type (DNARecord, RNARecord, or ProteinRecord)",
                ));
            }
            kind = Some(RecordKind::Dna);
            dna_records.push(record.inner.clone());
            continue;
        }
        if let Ok(record) = item.extract::<PyRef<'_, RNARecord>>() {
            if matches!(kind, Some(RecordKind::Dna | RecordKind::Protein)) {
                return Err(PyTypeError::new_err(
                    "write_fasta records must all be the same type (DNARecord, RNARecord, or ProteinRecord)",
                ));
            }
            kind = Some(RecordKind::Rna);
            rna_records.push(record.inner.clone());
            continue;
        }
        if let Ok(record) = item.extract::<PyRef<'_, ProteinRecord>>() {
            if matches!(kind, Some(RecordKind::Dna | RecordKind::Rna)) {
                return Err(PyTypeError::new_err(
                    "write_fasta records must all be the same type (DNARecord, RNARecord, or ProteinRecord)",
                ));
            }
            kind = Some(RecordKind::Protein);
            protein_records.push(record.inner.clone());
            continue;
        }
        return Err(PyTypeError::new_err(
            "write_fasta expects DNARecordBatch, RNARecordBatch, ProteinRecordBatch, or an iterable of DNARecord, RNARecord, or ProteinRecord",
        ));
    }

    if !dna_records.is_empty() {
        return fasta::write_fasta_records_to_path(path, &dna_records, line_width)
            .map_err(map_bio_err);
    }
    if !rna_records.is_empty() {
        return fasta::write_fasta_records_to_path(path, &rna_records, line_width)
            .map_err(map_bio_err);
    }
    if !protein_records.is_empty() {
        return fasta::write_fasta_records_to_path(path, &protein_records, line_width)
            .map_err(map_bio_err);
    }

    Err(PyTypeError::new_err(
        "write_fasta expects a non-empty iterable of DNARecord, RNARecord, or ProteinRecord",
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
