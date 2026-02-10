#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyIOError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::dna_record_batch::DNARecordBatch;
use crate::protein_record_batch::ProteinRecordBatch;
use crate::report::SkippedRecord;
use crate::rna_record_batch::RNARecordBatch;
use biorust_core::error::BioError;
use biorust_core::io::csv as core_csv;
use biorust_core::io::csv::ColumnSel;
use biorust_core::io::detect::{detect_seq_type, SeqType};
use biorust_core::io::OnError;

#[pyfunction]
#[pyo3(signature = (path, *, id_col, seq_col, desc_col=None, alphabet="auto", on_error="raise"))]
fn read_csv(
    py: Python<'_>,
    path: &str,
    id_col: &Bound<'_, PyAny>,
    seq_col: &Bound<'_, PyAny>,
    desc_col: Option<&Bound<'_, PyAny>>,
    alphabet: &str,
    on_error: &str,
) -> PyResult<PyObject> {
    let id_col = parse_col_sel(id_col)?;
    let seq_col_sel = parse_col_sel(seq_col)?;
    let desc_col = match desc_col {
        Some(obj) => Some(parse_col_sel(obj)?),
        None => None,
    };

    let on_error = parse_on_error(on_error)?;
    let alpha = match alphabet.to_ascii_lowercase().as_str() {
        "auto" => detect_csv_type(path, &seq_col_sel)?,
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
            let report = py
                .allow_threads(|| {
                    core_csv::read_csv_dna(&path, id_col, seq_col_sel, desc_col, on_error)
                })
                .map_err(map_bio_err)?;
            let skipped = report
                .skipped
                .into_iter()
                .map(SkippedRecord::from)
                .collect();
            let out = DNARecordBatch {
                inner: report.data,
                skipped,
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
        SeqType::Rna => {
            let report = py
                .allow_threads(|| {
                    core_csv::read_csv_rna(&path, id_col, seq_col_sel, desc_col, on_error)
                })
                .map_err(map_bio_err)?;
            let skipped = report
                .skipped
                .into_iter()
                .map(SkippedRecord::from)
                .collect();
            let out = RNARecordBatch {
                inner: report.data,
                skipped,
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
        SeqType::Protein => {
            let report = py
                .allow_threads(|| {
                    core_csv::read_csv_protein(&path, id_col, seq_col_sel, desc_col, on_error)
                })
                .map_err(map_bio_err)?;
            let skipped = report
                .skipped
                .into_iter()
                .map(SkippedRecord::from)
                .collect();
            let out = ProteinRecordBatch {
                inner: report.data,
                skipped,
            };
            Ok(Py::new(py, out)?.to_object(py))
        }
    }
}

/// Peek at the first ~100 rows of the seq column to detect the alphabet.
fn detect_csv_type(path: &str, seq_col: &ColumnSel) -> PyResult<SeqType> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .from_path(path)
        .map_err(|e| PyIOError::new_err(e.to_string()))?;

    let headers = rdr
        .headers()
        .map_err(|e| PyIOError::new_err(e.to_string()))?;
    if headers.is_empty() {
        return Err(PyValueError::new_err("CSV file is empty"));
    }

    let seq_idx = match seq_col {
        ColumnSel::Name(name) => headers
            .iter()
            .position(|h| h.trim() == name)
            .ok_or_else(|| PyValueError::new_err(format!("column '{name}' not found")))?,
        ColumnSel::Index(idx) => *idx,
    };

    let mut sample = Vec::new();
    for result in rdr.records().take(100) {
        let record = result.map_err(|e| PyIOError::new_err(e.to_string()))?;
        if let Some(field) = record.get(seq_idx) {
            for b in field.bytes() {
                if !b.is_ascii_whitespace() {
                    sample.push(b);
                }
            }
        }
    }

    Ok(detect_seq_type(&sample))
}

#[pyfunction]
fn csv_columns(path: &str) -> PyResult<Vec<String>> {
    core_csv::csv_columns(path).map_err(map_bio_err)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_csv, m)?)?;
    m.add_function(wrap_pyfunction!(csv_columns, m)?)?;
    Ok(())
}

fn parse_col_sel(obj: &Bound<'_, PyAny>) -> PyResult<ColumnSel> {
    if let Ok(name) = obj.extract::<String>() {
        return Ok(ColumnSel::Name(name));
    }
    if let Ok(index) = obj.extract::<isize>() {
        if index < 0 {
            return Err(PyValueError::new_err(
                "column index must be a non-negative integer",
            ));
        }
        return Ok(ColumnSel::Index(index as usize));
    }
    Err(PyTypeError::new_err(
        "column selector must be a string or non-negative integer",
    ))
}

fn map_bio_err(err: BioError) -> PyErr {
    match err {
        BioError::CsvParse { ref source, .. } if source.is_io_error() => {
            PyIOError::new_err(err.to_string())
        }
        BioError::FastaIo(_) => PyIOError::new_err(err.to_string()),
        _ => PyValueError::new_err(err.to_string()),
    }
}

fn parse_on_error(value: &str) -> PyResult<OnError> {
    match value.to_ascii_lowercase().as_str() {
        "raise" => Ok(OnError::Raise),
        "skip" => Ok(OnError::Skip),
        _ => Err(PyTypeError::new_err("on_error must be 'raise' or 'skip'")),
    }
}
