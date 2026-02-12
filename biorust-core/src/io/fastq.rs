use crate::error::{BioError, BioResult};
use crate::seq::record::SeqRecord;
use crate::seq::record_batch::RecordBatch;
use crate::seq::traits::SeqBytes;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Cursor, Write};
use std::marker::PhantomData;
use std::path::Path;

pub struct FastqRecords<R, S> {
    reader: R,
    line_no: usize,
    buf_line: String,
    _marker: PhantomData<S>,
}

impl<R: BufRead, S: SeqBytes> FastqRecords<R, S> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_no: 0,
            buf_line: String::new(),
            _marker: PhantomData,
        }
    }

    fn next_nonempty_line(&mut self) -> Option<BioResult<(String, usize)>> {
        loop {
            let (line, line_no) = match self.next_line() {
                Some(Ok(value)) => value,
                Some(Err(err)) => return Some(Err(err)),
                None => return None,
            };
            if line.trim().is_empty() {
                continue;
            }
            return Some(Ok((line, line_no)));
        }
    }

    fn next_line(&mut self) -> Option<BioResult<(String, usize)>> {
        self.buf_line.clear();
        match self.reader.read_line(&mut self.buf_line) {
            Ok(0) => None,
            Ok(_) => {
                self.line_no += 1;
                let line_no = self.line_no;
                Some(Ok((std::mem::take(&mut self.buf_line), line_no)))
            }
            Err(err) => Some(Err(BioError::FastqIo(err))),
        }
    }

    fn read_required_line(&mut self, msg: &'static str, line: usize) -> BioResult<(String, usize)> {
        match self.next_line() {
            Some(Ok(value)) => Ok(value),
            Some(Err(err)) => Err(err),
            None => Err(BioError::FastqFormat { msg, line }),
        }
    }
}

impl<R: BufRead, S: SeqBytes> Iterator for FastqRecords<R, S> {
    type Item = BioResult<SeqRecord<S>>;

    fn next(&mut self) -> Option<Self::Item> {
        let (header_line, header_line_no) = match self.next_nonempty_line()? {
            Ok(value) => value,
            Err(err) => return Some(Err(err)),
        };

        if !header_line.starts_with('@') {
            return Some(Err(BioError::FastqFormat {
                msg: "expected header line starting with '@'",
                line: header_line_no,
            }));
        }

        let (id, desc) = match parse_header(&header_line, header_line_no) {
            Ok(parsed) => parsed,
            Err(err) => return Some(Err(err)),
        };

        let (seq_line, seq_line_no) = match self
            .read_required_line("missing sequence line", header_line_no.saturating_add(1))
        {
            Ok(value) => value,
            Err(err) => return Some(Err(err)),
        };

        let (plus_line, plus_line_no) = match self
            .read_required_line("missing '+' separator line", seq_line_no.saturating_add(1))
        {
            Ok(value) => value,
            Err(err) => return Some(Err(err)),
        };

        if !plus_line.starts_with('+') {
            return Some(Err(BioError::FastqFormat {
                msg: "expected '+' separator line",
                line: plus_line_no,
            }));
        }

        let (qual_line, qual_line_no) =
            match self.read_required_line("missing quality line", plus_line_no.saturating_add(1)) {
                Ok(value) => value,
                Err(err) => return Some(Err(err)),
            };

        let seq_line = trim_eol(&seq_line);
        let qual_line = trim_eol(&qual_line);
        let seq_bytes = seq_line.as_bytes().to_vec();

        if seq_bytes.len() != qual_line.len() {
            return Some(Err(BioError::FastqFormat {
                msg: "sequence and quality lengths differ",
                line: qual_line_no,
            }));
        }

        let seq = match S::from_bytes(seq_bytes) {
            Ok(seq) => seq,
            Err(err) => return Some(Err(err)),
        };

        let record = match desc {
            Some(desc) => SeqRecord::new(id, seq).with_desc(desc),
            None => SeqRecord::new(id, seq),
        };
        Some(Ok(record))
    }
}

pub fn fastq_records_from_reader<R: BufRead, S: SeqBytes>(reader: R) -> FastqRecords<R, S> {
    FastqRecords::new(reader)
}

pub fn read_fastq_records_from_reader<R: BufRead, S: SeqBytes>(
    reader: R,
) -> BioResult<Vec<SeqRecord<S>>> {
    let mut out = Vec::new();
    for record in fastq_records_from_reader(reader) {
        out.push(record?);
    }
    Ok(out)
}

pub fn read_fastq_records_from_path<S: SeqBytes>(
    path: impl AsRef<Path>,
) -> BioResult<Vec<SeqRecord<S>>> {
    let file = File::open(path).map_err(BioError::FastqIo)?;
    let reader = BufReader::new(file);
    read_fastq_records_from_reader(reader)
}

pub fn read_fastq_records_from_bytes<S: SeqBytes>(data: &[u8]) -> BioResult<Vec<SeqRecord<S>>> {
    let reader = BufReader::new(Cursor::new(data));
    read_fastq_records_from_reader(reader)
}

pub fn read_fastq_batch_from_reader<R: BufRead, S: SeqBytes>(
    reader: R,
) -> BioResult<RecordBatch<S>> {
    let records = read_fastq_records_from_reader(reader)?;
    Ok(RecordBatch::from_records(records))
}

pub fn read_fastq_batch_from_path<S: SeqBytes>(
    path: impl AsRef<Path>,
) -> BioResult<RecordBatch<S>> {
    let file = File::open(path).map_err(BioError::FastqIo)?;
    let reader = BufReader::new(file);
    read_fastq_batch_from_reader(reader)
}

pub fn read_fastq_batch_from_bytes<S: SeqBytes>(data: &[u8]) -> BioResult<RecordBatch<S>> {
    let reader = BufReader::new(Cursor::new(data));
    read_fastq_batch_from_reader(reader)
}

pub fn write_fastq_records_to_writer<W: Write, S: SeqBytes>(
    writer: W,
    records: &[SeqRecord<S>],
    quality_char: u8,
) -> BioResult<()> {
    validate_quality_char(quality_char)?;
    let mut writer = BufWriter::new(writer);
    for record in records {
        write_fastq_record(
            &mut writer,
            &record.id,
            record.desc.as_deref(),
            record.seq.as_bytes(),
            quality_char,
        )?;
    }
    writer.flush().map_err(BioError::FastqIo)?;
    Ok(())
}

pub fn write_fastq_records_to_path<S: SeqBytes>(
    path: impl AsRef<Path>,
    records: &[SeqRecord<S>],
    quality_char: u8,
) -> BioResult<()> {
    let file = File::create(path).map_err(BioError::FastqIo)?;
    write_fastq_records_to_writer(file, records, quality_char)
}

pub fn write_fastq_batch_to_writer<W: Write, S: SeqBytes>(
    writer: W,
    batch: &RecordBatch<S>,
    quality_char: u8,
) -> BioResult<()> {
    validate_quality_char(quality_char)?;
    let mut writer = BufWriter::new(writer);
    for i in 0..batch.len() {
        let id = batch.id(i).expect("record batch length is consistent");
        let desc = batch.descs().get(i).and_then(|d| d.as_deref());
        let seq = batch
            .seq(i)
            .expect("record batch length is consistent")
            .as_bytes();
        write_fastq_record(&mut writer, id, desc, seq, quality_char)?;
    }
    writer.flush().map_err(BioError::FastqIo)?;
    Ok(())
}

pub fn write_fastq_batch_to_path<S: SeqBytes>(
    path: impl AsRef<Path>,
    batch: &RecordBatch<S>,
    quality_char: u8,
) -> BioResult<()> {
    let file = File::create(path).map_err(BioError::FastqIo)?;
    write_fastq_batch_to_writer(file, batch, quality_char)
}

fn validate_quality_char(ch: u8) -> BioResult<()> {
    if ch == b'\n' || ch == b'\r' {
        return Err(BioError::FastqInvalidQualityChar { ch: ch as char });
    }
    Ok(())
}

fn parse_header(header_line: &str, line_no: usize) -> BioResult<(Box<str>, Option<Box<str>>)> {
    let header = header_line.strip_prefix('@').ok_or(BioError::FastqFormat {
        msg: "expected header line starting with '@'",
        line: line_no,
    })?;

    let header = trim_eol(header).trim_start();
    if header.is_empty() {
        return Err(BioError::FastqFormat {
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

fn write_fastq_record<W: Write>(
    writer: &mut W,
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
    quality_char: u8,
) -> BioResult<()> {
    writer.write_all(b"@").map_err(BioError::FastqIo)?;
    write_header_field(writer, id)?;
    if let Some(desc) = desc {
        if !desc.is_empty() {
            writer.write_all(b" ").map_err(BioError::FastqIo)?;
            write_header_field(writer, desc)?;
        }
    }
    writer.write_all(b"\n").map_err(BioError::FastqIo)?;
    writer.write_all(seq).map_err(BioError::FastqIo)?;
    writer.write_all(b"\n+\n").map_err(BioError::FastqIo)?;
    let qual = vec![quality_char; seq.len()];
    writer.write_all(&qual).map_err(BioError::FastqIo)?;
    writer.write_all(b"\n").map_err(BioError::FastqIo)?;
    Ok(())
}

fn write_header_field<W: Write>(writer: &mut W, value: &str) -> BioResult<()> {
    let bytes = value.as_bytes();
    if bytes.contains(&b'\n') || bytes.contains(&b'\r') {
        for &b in bytes {
            match b {
                b'\n' | b'\r' => writer.write_all(b" ").map_err(BioError::FastqIo)?,
                _ => writer.write_all(&[b]).map_err(BioError::FastqIo)?,
            }
        }
    } else {
        writer.write_all(bytes).map_err(BioError::FastqIo)?;
    }
    Ok(())
}

fn trim_eol(line: &str) -> &str {
    line.trim_end_matches(&['\n', '\r'][..])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::DnaSeq;

    #[test]
    fn parse_single_record() {
        let data = b"@seq1\nACGT\n+\n!!!!\n";
        let records = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id(), "seq1");
        assert_eq!(records[0].desc(), None);
        assert_eq!(records[0].seq().as_bytes(), b"ACGT");
    }

    #[test]
    fn header_with_description() {
        let data = b"@seq1 some desc here\nACGT\n+\nIIII\n";
        let records = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records[0].id(), "seq1");
        assert_eq!(records[0].desc(), Some("some desc here"));
    }

    #[test]
    fn multiple_records() {
        let data = b"@seq1\nAC\n+\n!!\n@seq2\nGT\n+\n##\n";
        let records = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), "seq1");
        assert_eq!(records[1].id(), "seq2");
    }

    #[test]
    fn empty_sequence_allowed() {
        let data = b"@seq1\n\n+\n\n";
        let records = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].seq().as_bytes(), b"");
    }

    #[test]
    fn invalid_header() {
        let data = b">seq1\nAC\n+\n!!\n";
        let err = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap_err();
        match err {
            BioError::FastqFormat { .. } => {}
            other => panic!("expected fastq format error, got {other:?}"),
        }
    }

    #[test]
    fn invalid_plus_separator() {
        let data = b"@seq1\nAC\n-\n!!\n";
        let err = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap_err();
        match err {
            BioError::FastqFormat { .. } => {}
            other => panic!("expected fastq format error, got {other:?}"),
        }
    }

    #[test]
    fn quality_length_mismatch() {
        let data = b"@seq1\nACGT\n+\n!!!\n";
        let err = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap_err();
        match err {
            BioError::FastqFormat { .. } => {}
            other => panic!("expected fastq format error, got {other:?}"),
        }
    }

    #[test]
    fn invalid_sequence_char() {
        let data = b"@seq1\nAC#\n+\n!!!\n";
        let err = read_fastq_records_from_bytes::<DnaSeq>(data).unwrap_err();
        match err {
            BioError::InvalidChar { .. } => {}
            other => panic!("expected invalid char error, got {other:?}"),
        }
    }

    #[test]
    fn write_records() {
        let records = vec![SeqRecord::new(
            "seq1",
            DnaSeq::new(b"ACGT".to_vec()).unwrap(),
        )];
        let mut out = Vec::new();
        write_fastq_records_to_writer(&mut out, &records, b'I').unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text, "@seq1\nACGT\n+\nIIII\n");
    }
}
