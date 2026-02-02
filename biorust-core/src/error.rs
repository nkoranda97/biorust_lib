use thiserror::Error;

#[derive(Debug, Error)]
pub enum BioError {
    #[error("invalid character '{ch}' at position {pos}")]
    InvalidChar { ch: char, pos: usize },

    #[error("invalid frame: {frame} (must be 0, 1, or 2)")]
    InvalidFrame { frame: usize },
}

pub type BioResult<T> = Result<T, BioError>;