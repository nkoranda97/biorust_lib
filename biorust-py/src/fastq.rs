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
use biorust_core::io::fastq as core_fastq;
use biorust_core::seq::dna::DnaSeq;
use biorust_core::seq::protein::ProteinSeq;
use biorust_core::seq::record::SeqRecord;
use biorust_core::seq::rna::RnaSeq;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[pyfunction]
#[pyo3(signature = (path, *, alphabet="auto"))]
fn read_fastq(py: Python<'_>, path: &str, alphabet: &str) -> PyResult<PyObject> {
    let alpha = match alphabet.to_ascii_lowercase().as_str() {
        "auto" => detect_fastq_type(path)?,
        "dna" => SeqType::Dna,
        "rna" => SeqType::Rna,
        "protein" => SeqType::Protein,
        _ => {
            return Err(PyValueError::new_err(
                "alphabet must be 'auto', 'dna', 'rna', or 'protein'",
            ))
        }
    };

    let path = path.to_owned();
    match alpha {
        SeqType::Dna => {
            let batch = py
                .allow_threads(|| core_fastq::read_fastq_batch_from_path::<DnaSeq>(&path))
                .map_err(map_bio_err)?;
            let out = DNARecordBatch {
                inner: batch,
                skipped: Vec::new(),
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
        SeqType::Rna => {
            let batch = py
                .allow_threads(|| core_fastq::read_fastq_batch_from_path::<RnaSeq>(&path))
                .map_err(map_bio_err)?;
            let out = RNARecordBatch {
                inner: batch,
                skipped: Vec::new(),
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
        SeqType::Protein => {
            let batch = py
                .allow_threads(|| core_fastq::read_fastq_batch_from_path::<ProteinSeq>(&path))
                .map_err(map_bio_err)?;
            let out = ProteinRecordBatch {
                inner: batch,
                skipped: Vec::new(),
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
    }
}

/// Peek at the first FASTQ record's sequence bytes to detect the alphabet.
fn detect_fastq_type(path: &str) -> PyResult<SeqType> {
    let file = File::open(path).map_err(|e| PyIOError::new_err(e.to_string()))?;
    let mut lines = BufReader::new(file).lines();

    loop {
        let line = match lines.next() {
            Some(Ok(line)) => line,
            Some(Err(err)) => return Err(PyIOError::new_err(err.to_string())),
            None => return Ok(SeqType::Dna),
        };

        if line.trim().is_empty() {
            continue;
        }
        if !line.starts_with('@') {
            return Err(PyValueError::new_err(
                "invalid FASTQ: expected header line starting with '@'",
            ));
        }

        let seq_line = match lines.next() {
            Some(Ok(line)) => line,
            Some(Err(err)) => return Err(PyIOError::new_err(err.to_string())),
            None => {
                return Err(PyValueError::new_err(
                    "invalid FASTQ: missing sequence line",
                ))
            }
        };
        let sample = seq_line
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .collect::<Vec<u8>>();
        return Ok(detect_seq_type(&sample));
    }
}

#[pyfunction]
#[pyo3(signature = (path, records, *, quality_char="I"))]
fn write_fastq(path: &str, records: &Bound<'_, PyAny>, quality_char: &str) -> PyResult<()> {
    let quality_char = parse_quality_char(quality_char)?;

    if let Ok(batch) = records.extract::<PyRef<'_, DNARecordBatch>>() {
        return core_fastq::write_fastq_batch_to_path(path, &batch.inner, quality_char)
            .map_err(map_bio_err);
    }
    if let Ok(batch) = records.extract::<PyRef<'_, ProteinRecordBatch>>() {
        return core_fastq::write_fastq_batch_to_path(path, &batch.inner, quality_char)
            .map_err(map_bio_err);
    }
    if let Ok(batch) = records.extract::<PyRef<'_, RNARecordBatch>>() {
        return core_fastq::write_fastq_batch_to_path(path, &batch.inner, quality_char)
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
                    "write_fastq records must all be the same type (DNARecord, RNARecord, or ProteinRecord)",
                ));
            }
            kind = Some(RecordKind::Dna);
            dna_records.push(record.inner.clone());
            continue;
        }
        if let Ok(record) = item.extract::<PyRef<'_, RNARecord>>() {
            if matches!(kind, Some(RecordKind::Dna | RecordKind::Protein)) {
                return Err(PyTypeError::new_err(
                    "write_fastq records must all be the same type (DNARecord, RNARecord, or ProteinRecord)",
                ));
            }
            kind = Some(RecordKind::Rna);
            rna_records.push(record.inner.clone());
            continue;
        }
        if let Ok(record) = item.extract::<PyRef<'_, ProteinRecord>>() {
            if matches!(kind, Some(RecordKind::Dna | RecordKind::Rna)) {
                return Err(PyTypeError::new_err(
                    "write_fastq records must all be the same type (DNARecord, RNARecord, or ProteinRecord)",
                ));
            }
            kind = Some(RecordKind::Protein);
            protein_records.push(record.inner.clone());
            continue;
        }
        return Err(PyTypeError::new_err(
            "write_fastq expects DNARecordBatch, RNARecordBatch, ProteinRecordBatch, or an iterable of DNARecord, RNARecord, or ProteinRecord",
        ));
    }

    if !dna_records.is_empty() {
        return core_fastq::write_fastq_records_to_path(path, &dna_records, quality_char)
            .map_err(map_bio_err);
    }
    if !rna_records.is_empty() {
        return core_fastq::write_fastq_records_to_path(path, &rna_records, quality_char)
            .map_err(map_bio_err);
    }
    if !protein_records.is_empty() {
        return core_fastq::write_fastq_records_to_path(path, &protein_records, quality_char)
            .map_err(map_bio_err);
    }

    Err(PyTypeError::new_err(
        "write_fastq expects a non-empty iterable of DNARecord, RNARecord, or ProteinRecord",
    ))
}

fn parse_quality_char(value: &str) -> PyResult<u8> {
    let mut chars = value.chars();
    let ch = chars
        .next()
        .ok_or_else(|| PyValueError::new_err("quality_char must be a single character"))?;
    if chars.next().is_some() {
        return Err(PyValueError::new_err(
            "quality_char must be a single character",
        ));
    }
    if ch == '\n' || ch == '\r' {
        return Err(PyValueError::new_err(
            "quality_char must not be a newline character",
        ));
    }
    if !ch.is_ascii() {
        return Err(PyValueError::new_err("quality_char must be ASCII"));
    }
    Ok(ch as u8)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(write_fastq, m)?)?;
    Ok(())
}

fn map_bio_err(err: BioError) -> PyErr {
    match err {
        BioError::FastqIo(io) => PyIOError::new_err(io.to_string()),
        other => PyValueError::new_err(other.to_string()),
    }
}
