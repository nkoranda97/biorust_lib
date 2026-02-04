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

    #[error("fasta format error at line {line}: {msg}")]
    FastaFormat { msg: &'static str, line: usize },

    #[error("fasta io error: {0}")]
    FastaIo(#[from] io::Error),
}

pub type BioResult<T> = Result<T, BioError>;
