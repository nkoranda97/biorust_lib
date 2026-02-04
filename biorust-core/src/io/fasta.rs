use crate::error::{BioError, BioResult};
use crate::seq::record::SeqRecord;
use crate::seq::record_batch::RecordBatch;
use crate::seq::traits::SeqBytes;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor};
use std::marker::PhantomData;
use std::path::Path;

pub struct FastaRecords<R, S> {
    reader: R,
    line_no: usize,
    pending_header: Option<(String, usize)>,
    buf_line: String,
    seq_buf: Vec<u8>,
    _marker: PhantomData<S>,
}

impl<R: BufRead, S: SeqBytes> FastaRecords<R, S> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_no: 0,
            pending_header: None,
            buf_line: String::new(),
            seq_buf: Vec::new(),
            _marker: PhantomData,
        }
    }

    fn next_header(&mut self) -> Option<BioResult<(String, usize)>> {
        if let Some(pending) = self.pending_header.take() {
            return Some(Ok(pending));
        }

        loop {
            self.buf_line.clear();
            match self.reader.read_line(&mut self.buf_line) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line_no += 1;
                    let line_no = self.line_no;
                    if self.buf_line.starts_with('>') {
                        return Some(Ok((self.buf_line.clone(), line_no)));
                    }
                    if self.buf_line.trim().is_empty() {
                        continue;
                    }
                    return Some(Err(BioError::FastaFormat {
                        msg: "expected header line starting with '>'",
                        line: line_no,
                    }));
                }
                Err(err) => return Some(Err(BioError::FastaIo(err))),
            }
        }
    }
}

impl<R: BufRead, S: SeqBytes> Iterator for FastaRecords<R, S> {
    type Item = BioResult<SeqRecord<S>>;

    fn next(&mut self) -> Option<Self::Item> {
        let (header_line, header_line_no) = match self.next_header()? {
            Ok(header) => header,
            Err(err) => return Some(Err(err)),
        };

        let (id, desc) = match parse_header(&header_line, header_line_no) {
            Ok(parsed) => parsed,
            Err(err) => return Some(Err(err)),
        };

        self.seq_buf.clear();

        loop {
            self.buf_line.clear();
            match self.reader.read_line(&mut self.buf_line) {
                Ok(0) => break,
                Ok(_) => {
                    self.line_no += 1;
                    let line_no = self.line_no;
                    if self.buf_line.starts_with('>') {
                        self.pending_header = Some((self.buf_line.clone(), line_no));
                        break;
                    }
                    for b in self.buf_line.bytes() {
                        if !b.is_ascii_whitespace() {
                            self.seq_buf.push(b);
                        }
                    }
                }
                Err(err) => return Some(Err(BioError::FastaIo(err))),
            }
        }

        let capacity = self.seq_buf.capacity();
        let bytes = std::mem::take(&mut self.seq_buf);
        let seq = match S::from_bytes(bytes) {
            Ok(seq) => seq,
            Err(err) => return Some(Err(err)),
        };
        self.seq_buf = Vec::with_capacity(capacity);

        Some(Ok(SeqRecord { id, desc, seq }))
    }
}

pub fn fasta_records_from_reader<R: BufRead, S: SeqBytes>(reader: R) -> FastaRecords<R, S> {
    FastaRecords::new(reader)
}

pub fn read_fasta_records_from_reader<R: BufRead, S: SeqBytes>(
    reader: R,
) -> BioResult<Vec<SeqRecord<S>>> {
    let mut out = Vec::new();
    for record in fasta_records_from_reader(reader) {
        out.push(record?);
    }
    Ok(out)
}

pub fn read_fasta_records_from_path<S: SeqBytes>(
    path: impl AsRef<Path>,
) -> BioResult<Vec<SeqRecord<S>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    read_fasta_records_from_reader(reader)
}

pub fn read_fasta_records_from_bytes<S: SeqBytes>(data: &[u8]) -> BioResult<Vec<SeqRecord<S>>> {
    let reader = BufReader::new(Cursor::new(data));
    read_fasta_records_from_reader(reader)
}

pub fn read_fasta_batch_from_reader<R: BufRead, S: SeqBytes>(
    reader: R,
) -> BioResult<RecordBatch<S>> {
    let records = read_fasta_records_from_reader(reader)?;
    Ok(RecordBatch::from_records(records))
}

pub fn read_fasta_batch_from_path<S: SeqBytes>(
    path: impl AsRef<Path>,
) -> BioResult<RecordBatch<S>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    read_fasta_batch_from_reader(reader)
}

pub fn read_fasta_batch_from_bytes<S: SeqBytes>(data: &[u8]) -> BioResult<RecordBatch<S>> {
    let reader = BufReader::new(Cursor::new(data));
    read_fasta_batch_from_reader(reader)
}

fn parse_header(header_line: &str, line_no: usize) -> BioResult<(Box<str>, Option<Box<str>>)> {
    let header = header_line.strip_prefix('>').ok_or(BioError::FastaFormat {
        msg: "expected header line starting with '>'",
        line: line_no,
    })?;

    let header = header.trim_end_matches(&['\n', '\r'][..]).trim_start();
    if header.is_empty() {
        return Err(BioError::FastaFormat {
            msg: "empty header",
            line: line_no,
        });
    }

    let (id, desc) = match header.find(|c: char| c.is_whitespace()) {
        Some(idx) => {
            let id = &header[..idx];
            let desc = header[idx..].trim();
            let desc = if desc.is_empty() { None } else { Some(desc) };
            (id, desc)
        }
        None => (header, None),
    };

    Ok((id.into(), desc.map(|s| s.into())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::DnaSeq;

    #[test]
    fn parse_single_record() {
        let data = b">seq1\nACGT\n";
        let records = read_fasta_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id(), "seq1");
        assert_eq!(records[0].desc(), None);
        assert_eq!(records[0].seq().as_bytes(), b"ACGT");
    }

    #[test]
    fn header_with_description() {
        let data = b">seq1 some desc here\nAC\nGT\n";
        let records = read_fasta_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records[0].id(), "seq1");
        assert_eq!(records[0].desc(), Some("some desc here"));
        assert_eq!(records[0].seq().as_bytes(), b"ACGT");
    }

    #[test]
    fn multiple_records() {
        let data = b">seq1\nAC\n>seq2\nGT\n";
        let records = read_fasta_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "seq1");
        assert_eq!(records[1].id(), "seq2");
    }

    #[test]
    fn empty_sequence_allowed() {
        let data = b">seq1\n>seq2\nA\n";
        let records = read_fasta_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].seq().as_bytes(), b"");
        assert_eq!(records[1].seq().as_bytes(), b"A");
    }

    #[test]
    fn invalid_format_before_header() {
        let data = b"ACGT\n>seq1\nAC\n";
        let err = read_fasta_records_from_bytes::<DnaSeq>(data).unwrap_err();
        match err {
            BioError::FastaFormat { .. } => {}
            other => panic!("expected fasta format error, got {other:?}"),
        }
    }

    #[test]
    fn invalid_sequence_char() {
        let data = b">seq1\nAC#\n";
        let err = read_fasta_records_from_bytes::<DnaSeq>(data).unwrap_err();
        match err {
            BioError::InvalidChar { .. } => {}
            other => panic!("expected invalid char error, got {other:?}"),
        }
    }
}
