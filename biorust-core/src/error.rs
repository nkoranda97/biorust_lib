use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum BioError {
    #[error("invalid character '{ch}' at position {pos}")]
    InvalidChar { ch: char, pos: usize },

    #[error("invalid frame: {frame} (must be 0, 1, or 2)")]
    InvalidFrame { frame: usize },

    #[error("integer byte out of range: {val} (expected 0..=255)")]
    IntByteOutOfRange { val: i128 },

    #[error("invalid window size: {window}")]
    InvalidWindow { window: usize },

    #[error("invalid scoring parameters: {msg}")]
    InvalidScoring { msg: String },

    #[error("fasta format error at line {line}: {msg}")]
    FastaFormat { msg: &'static str, line: usize },

    #[error("fasta io error: {0}")]
    FastaIo(#[from] io::Error),

    #[error("record batch length mismatch (ids={ids}, descs={descs}, seqs={seqs})")]
    RecordBatchLenMismatch {
        ids: usize,
        descs: usize,
        seqs: usize,
    },

    #[error("csv missing column '{name}' in {path}. headers: {headers:?}")]
    CsvMissingColumn {
        name: String,
        headers: Vec<String>,
        path: String,
    },

    #[error("csv column index {index} out of range (ncols={ncols}) in {path}")]
    CsvColumnIndexOutOfRange {
        index: usize,
        ncols: usize,
        path: String,
    },

    #[error("csv missing field at row {row} for column {column} in {path}")]
    CsvMissingField {
        row: usize,
        column: String,
        path: String,
    },

    #[error("csv invalid sequence at row {row} for column {column} in {path}: {source}")]
    CsvInvalidSequence {
        row: usize,
        column: String,
        path: String,
        #[source]
        source: Box<BioError>,
    },

    #[error("csv parse error in {path}: {source}")]
    CsvParse {
        path: String,
        #[source]
        source: csv::Error,
    },
}

pub type BioResult<T> = Result<T, BioError>;
