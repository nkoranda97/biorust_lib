use crate::error::{BioError, BioResult};
use crate::io::{normalize_seq_bytes, OnError, ReadReport, SkippedRecord};
use crate::seq::dna::DnaSeq;
use crate::seq::protein::ProteinSeq;
use crate::seq::record_batch::RecordBatch;
use crate::seq::traits::SeqBytes;
use csv::{ReaderBuilder, StringRecord};
use std::fs::File;
use std::path::Path;

#[derive(Clone, Debug)]
pub enum ColumnSel {
    Name(String),
    Index(usize),
}

pub fn csv_columns(path: impl AsRef<Path>) -> BioResult<Vec<String>> {
    let path_ref = path.as_ref();
    let path_str = path_ref.display().to_string();
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path_ref)
        .map_err(|e| BioError::CsvParse {
            path: path_str.clone(),
            source: e,
        })?;
    let headers = reader
        .headers()
        .map_err(|e| BioError::CsvParse {
            path: path_str.clone(),
            source: e,
        })?
        .clone();
    Ok(headers.iter().map(|s| s.to_string()).collect())
}

pub fn read_csv_dna(
    path: impl AsRef<Path>,
    id_col: ColumnSel,
    seq_col: ColumnSel,
    desc_col: Option<ColumnSel>,
    on_error: OnError,
) -> BioResult<ReadReport<RecordBatch<DnaSeq>>> {
    read_csv(path, id_col, seq_col, desc_col, on_error)
}

pub fn read_csv_protein(
    path: impl AsRef<Path>,
    id_col: ColumnSel,
    seq_col: ColumnSel,
    desc_col: Option<ColumnSel>,
    on_error: OnError,
) -> BioResult<ReadReport<RecordBatch<ProteinSeq>>> {
    read_csv(path, id_col, seq_col, desc_col, on_error)
}

pub fn read_csv<S: SeqBytes>(
    path: impl AsRef<Path>,
    id_col: ColumnSel,
    seq_col: ColumnSel,
    desc_col: Option<ColumnSel>,
    on_error: OnError,
) -> BioResult<ReadReport<RecordBatch<S>>> {
    let path_ref = path.as_ref();
    let path_str = path_ref.display().to_string();
    let file = File::open(path_ref).map_err(|e| BioError::CsvParse {
        path: path_str.clone(),
        source: csv::Error::from(e),
    })?;

    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_reader(file);

    let headers = reader
        .headers()
        .map_err(|e| BioError::CsvParse {
            path: path_str.clone(),
            source: e,
        })?
        .clone();
    let id_idx = resolve_column(&id_col, &headers, &path_str)?;
    let seq_idx = resolve_column(&seq_col, &headers, &path_str)?;
    let desc_idx = desc_col
        .as_ref()
        .map(|sel| resolve_column(sel, &headers, &path_str))
        .transpose()?;

    let mut ids: Vec<Box<str>> = Vec::new();
    let mut descs: Vec<Option<Box<str>>> = Vec::new();
    let mut seqs: Vec<S> = Vec::new();
    let mut skipped: Vec<SkippedRecord> = Vec::new();

    for (row_idx, result) in reader.records().enumerate() {
        let record = result.map_err(|e| BioError::CsvParse {
            path: path_str.clone(),
            source: e,
        })?;
        let row = row_idx + 1;

        let id_field = record
            .get(id_idx)
            .ok_or_else(|| BioError::CsvMissingField {
                row,
                column: column_label(&id_col),
                path: path_str.clone(),
            })?;
        let id_value = id_field.trim();

        let seq_field = record
            .get(seq_idx)
            .ok_or_else(|| BioError::CsvMissingField {
                row,
                column: column_label(&seq_col),
                path: path_str.clone(),
            })?;
        let seq_bytes = normalize_seq_bytes(seq_field);
        match S::from_bytes(seq_bytes) {
            Ok(seq) => {
                ids.push(id_value.to_string().into_boxed_str());
                seqs.push(seq);
            }
            Err(err) => match on_error {
                OnError::Raise => {
                    return Err(BioError::CsvInvalidSequence {
                        row,
                        column: column_label(&seq_col),
                        path: path_str.clone(),
                        source: Box::new(err),
                    });
                }
                OnError::Skip => {
                    let msg = format!(
                        "invalid sequence at row {row}, column {}: {err}",
                        column_label(&seq_col)
                    );
                    let id = if id_value.is_empty() {
                        None
                    } else {
                        Some(id_value.to_string().into_boxed_str())
                    };
                    skipped.push(SkippedRecord {
                        row,
                        id,
                        column: column_label(&seq_col).into_boxed_str(),
                        message: msg.into_boxed_str(),
                    });
                    continue;
                }
            },
        }

        if let Some(desc_idx) = desc_idx {
            let desc_field = record
                .get(desc_idx)
                .ok_or_else(|| BioError::CsvMissingField {
                    row,
                    column: column_label(desc_col.as_ref().expect("desc_idx exists")),
                    path: path_str.clone(),
                })?;
            let desc = desc_field.trim();
            if desc.is_empty() {
                descs.push(None);
            } else {
                descs.push(Some(desc.to_string().into_boxed_str()));
            }
        } else {
            descs.push(None);
        }
    }

    let data = RecordBatch::new(ids, descs, seqs)?;
    Ok(ReadReport { data, skipped })
}

fn resolve_column(sel: &ColumnSel, headers: &StringRecord, path: &str) -> BioResult<usize> {
    match sel {
        ColumnSel::Name(name) => {
            headers
                .iter()
                .position(|h| h == name)
                .ok_or_else(|| BioError::CsvMissingColumn {
                    name: name.clone(),
                    headers: headers.iter().map(|s| s.to_string()).collect(),
                    path: path.to_string(),
                })
        }
        ColumnSel::Index(index) => {
            if *index < headers.len() {
                Ok(*index)
            } else {
                Err(BioError::CsvColumnIndexOutOfRange {
                    index: *index,
                    ncols: headers.len(),
                    path: path.to_string(),
                })
            }
        }
    }
}

fn column_label(sel: &ColumnSel) -> String {
    match sel {
        ColumnSel::Name(name) => name.clone(),
        ColumnSel::Index(index) => format!("#{index}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::DnaSeq;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn write_temp_csv(contents: &str) -> std::path::PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!("biorust_csv_test_{nanos}.csv"));
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn read_csv_basic() {
        let path = write_temp_csv("id,seq,desc\ns1,ACGT,first\ns2,TT,");
        let report = read_csv::<DnaSeq>(
            &path,
            ColumnSel::Name("id".to_string()),
            ColumnSel::Name("seq".to_string()),
            Some(ColumnSel::Name("desc".to_string())),
            OnError::Raise,
        )
        .unwrap();
        assert_eq!(report.data.len(), 2);
        assert_eq!(report.data.ids()[0].as_ref(), "s1");
        assert_eq!(report.data.descs()[0].as_deref(), Some("first"));
        assert_eq!(report.data.seqs().as_slice()[0].as_bytes(), b"ACGT");
        assert!(report.skipped.is_empty());
    }

    #[test]
    fn missing_column_name() {
        let path = write_temp_csv("id,seq\ns1,ACGT\n");
        let err = read_csv::<DnaSeq>(
            &path,
            ColumnSel::Name("id".to_string()),
            ColumnSel::Name("missing".to_string()),
            None,
            OnError::Raise,
        )
        .unwrap_err();
        match err {
            BioError::CsvMissingColumn { .. } => {}
            other => panic!("expected missing column error, got {other:?}"),
        }
    }

    #[test]
    fn column_index_out_of_range() {
        let path = write_temp_csv("id,seq\ns1,ACGT\n");
        let err = read_csv::<DnaSeq>(
            &path,
            ColumnSel::Index(0),
            ColumnSel::Index(5),
            None,
            OnError::Raise,
        )
        .unwrap_err();
        match err {
            BioError::CsvColumnIndexOutOfRange { .. } => {}
            other => panic!("expected column index error, got {other:?}"),
        }
    }

    #[test]
    fn invalid_sequence_char() {
        let path = write_temp_csv("id,seq\ns1,AC#\n");
        let err = read_csv::<DnaSeq>(
            &path,
            ColumnSel::Name("id".to_string()),
            ColumnSel::Name("seq".to_string()),
            None,
            OnError::Raise,
        )
        .unwrap_err();
        match err {
            BioError::CsvInvalidSequence { .. } => {}
            other => panic!("expected invalid sequence error, got {other:?}"),
        }
    }

    #[test]
    fn missing_field_error() {
        let path = write_temp_csv("id,seq,desc\ns1,ACGT\n");
        let err = read_csv::<DnaSeq>(
            &path,
            ColumnSel::Name("id".to_string()),
            ColumnSel::Name("seq".to_string()),
            Some(ColumnSel::Name("desc".to_string())),
            OnError::Raise,
        )
        .unwrap_err();
        match err {
            BioError::CsvMissingField { .. } => {}
            other => panic!("expected missing field error, got {other:?}"),
        }
    }

    #[test]
    fn skip_invalid_sequence() {
        let path = write_temp_csv("id,seq\ns1,ACGT\ns2,AC#\ns3,TT\n");
        let report = read_csv::<DnaSeq>(
            &path,
            ColumnSel::Name("id".to_string()),
            ColumnSel::Name("seq".to_string()),
            None,
            OnError::Skip,
        )
        .unwrap();
        assert_eq!(report.data.len(), 2);
        assert_eq!(report.skipped.len(), 1);
        assert_eq!(report.skipped[0].row, 2);
        assert_eq!(report.skipped[0].id.as_deref(), Some("s2"));
        assert!(report.skipped[0].message.contains("invalid sequence"));
    }
}
