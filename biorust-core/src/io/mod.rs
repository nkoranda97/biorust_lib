pub mod csv;
pub mod detect;
pub mod fasta;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OnError {
    Raise,
    Skip,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SkippedRecord {
    pub row: usize,
    pub id: Option<Box<str>>,
    pub column: Box<str>,
    pub message: Box<str>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadReport<T> {
    pub data: T,
    pub skipped: Vec<SkippedRecord>,
}

pub fn normalize_seq_bytes(input: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(input.len());
    for b in input.bytes() {
        if !b.is_ascii_whitespace() {
            out.push(b);
        }
    }
    out
}
