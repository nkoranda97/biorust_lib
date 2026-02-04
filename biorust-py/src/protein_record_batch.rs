#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule, PySlice};

use crate::batch::ProteinBatch;
use crate::protein_record::ProteinRecord;
use crate::report::SkippedRecord;
use biorust_core::seq::batch::SeqBatch;
use biorust_core::seq::protein::ProteinSeq;
use biorust_core::seq::record::SeqRecord;
use biorust_core::seq::record_batch::RecordBatch;

#[pyclass]
pub struct ProteinRecordBatch {
    pub(crate) inner: RecordBatch<ProteinSeq>,
    pub(crate) skipped: Vec<SkippedRecord>,
}

fn collect_records(obj: &Bound<'_, PyAny>) -> PyResult<Vec<SeqRecord<ProteinSeq>>> {
    if let Ok(batch) = obj.extract::<PyRef<'_, ProteinRecordBatch>>() {
        let ids = batch.inner.ids().to_vec();
        let descs = batch.inner.descs().to_vec();
        let seqs = batch.inner.seqs().as_slice().to_vec();
        let mut out = Vec::with_capacity(seqs.len());
        for i in 0..seqs.len() {
            out.push(SeqRecord {
                id: ids[i].clone(),
                desc: descs[i].clone(),
                seq: seqs[i].clone(),
            });
        }
        return Ok(out);
    }

    let mut out = Vec::new();
    for item in obj.iter()? {
        let item = item?;
        let record = item.extract::<PyRef<'_, ProteinRecord>>().map_err(|_| {
            PyTypeError::new_err("ProteinRecordBatch expects ProteinRecord objects only")
        })?;
        out.push(record.inner.clone());
    }
    Ok(out)
}

#[pymethods]
impl ProteinRecordBatch {
    #[new]
    fn new(records: &Bound<'_, PyAny>) -> PyResult<Self> {
        let records = collect_records(records)?;
        Ok(Self {
            inner: RecordBatch::from_records(records),
            skipped: Vec::new(),
        })
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        if let Ok(slice) = index.downcast::<PySlice>() {
            let idx = slice.indices(self.inner.len() as isize)?;
            let (start, stop, step) = (idx.start, idx.stop, idx.step);
            let mut ids = Vec::new();
            let mut descs = Vec::new();
            let mut seqs = Vec::new();

            if step > 0 {
                let mut i = start;
                while i < stop {
                    let idx = i as usize;
                    ids.push(self.inner.ids()[idx].clone());
                    descs.push(self.inner.descs()[idx].clone());
                    seqs.push(self.inner.seqs().as_slice()[idx].clone());
                    i += step;
                }
            } else {
                let mut i = start;
                while i > stop {
                    let idx = i as usize;
                    ids.push(self.inner.ids()[idx].clone());
                    descs.push(self.inner.descs()[idx].clone());
                    seqs.push(self.inner.seqs().as_slice()[idx].clone());
                    i += step;
                }
            }

            let batch = ProteinRecordBatch {
                inner: RecordBatch::new(ids, descs, seqs)
                    .map_err(|e| PyTypeError::new_err(e.to_string()))?,
                skipped: self.skipped.clone(),
            };
            return Ok(Py::new(py, batch)?.to_object(py));
        }

        let index: isize = index
            .extract()
            .map_err(|_| PyTypeError::new_err("index must be int or slice"))?;
        let n = self.inner.len() as isize;
        let i = if index < 0 { index + n } else { index };

        if i < 0 || i >= n {
            return Err(PyIndexError::new_err("index out of range"));
        }

        let i = i as usize;
        let record = ProteinRecord {
            inner: SeqRecord {
                id: self.inner.ids()[i].clone(),
                desc: self.inner.descs()[i].clone(),
                seq: self.inner.seqs().as_slice()[i].clone(),
            },
        };
        Ok(Py::new(py, record)?.to_object(py))
    }

    fn ids(&self) -> Vec<String> {
        self.inner.ids().iter().map(|s| s.to_string()).collect()
    }

    #[getter]
    fn skipped(&self) -> Vec<SkippedRecord> {
        self.skipped.clone()
    }

    fn descriptions(&self) -> Vec<Option<String>> {
        self.inner
            .descs()
            .iter()
            .map(|d| d.as_deref().map(|s| s.to_string()))
            .collect()
    }

    fn seqs(&self) -> ProteinBatch {
        let seqs: Vec<ProteinSeq> = self.inner.seqs().as_slice().to_vec();
        ProteinBatch {
            inner: SeqBatch::new(seqs),
        }
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ProteinRecordBatch>()?;
    Ok(())
}
