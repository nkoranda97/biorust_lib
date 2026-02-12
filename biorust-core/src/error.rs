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

    #[error("invalid feature location: start {start} > end {end}")]
    InvalidLocation { start: usize, end: usize },

    #[error("invalid strand: {strand} (expected -1 or 1)")]
    InvalidStrand { strand: i8 },

    #[error("invalid feature type: empty")]
    InvalidFeatureType,

    #[error("fasta format error at line {line}: {msg}")]
    FastaFormat { msg: &'static str, line: usize },

    #[error("fasta io error: {0}")]
    FastaIo(#[from] io::Error),

    #[error("fastq format error at line {line}: {msg}")]
    FastqFormat { msg: &'static str, line: usize },

    #[error("fastq io error: {0}")]
    FastqIo(io::Error),

    #[error("invalid fastq quality character: {ch:?}")]
    FastqInvalidQualityChar { ch: char },

    #[error("record batch length mismatch (ids={ids}, descs={descs}, seqs={seqs})")]
    RecordBatchLenMismatch {
        ids: usize,
        descs: usize,
        seqs: usize,
    },

    #[error("batch index {index} out of range (len={len})")]
    BatchIndexOutOfRange { index: usize, len: usize },

    #[error("empty batch")]
    EmptyBatch,

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

    #[error("too few sequences: {n} (need at least 2)")]
    TooFewSequences { n: usize },

    #[error("saturated distance between sequences {i} and {j} for model {model}")]
    SaturatedDistance { i: usize, j: usize, model: String },

    #[error("no valid sites between sequences {i} and {j}")]
    NoValidSites { i: usize, j: usize },

    #[error("label count {labels} does not match sequence count {seqs}")]
    LabelCountMismatch { labels: usize, seqs: usize },

    #[error("sequence {index} has length {len} but expected {expected}")]
    SequenceLengthMismatch {
        index: usize,
        len: usize,
        expected: usize,
    },

    #[error("translation error: {msg}")]
    TranslationError { msg: String },
}

pub type BioResult<T> = Result<T, BioError>;
