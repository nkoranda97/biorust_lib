use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::protein::Protein;
use crate::seq_shared;
use biorust_core::seq::protein::ProteinSeq;
use biorust_core::seq::record::SeqRecord;

#[pyclass(frozen)]
pub struct ProteinRecord {
    pub(crate) inner: SeqRecord<ProteinSeq>,
}

#[pymethods]
impl ProteinRecord {
    #[new]
    #[pyo3(signature = (id, seq, desc=None))]
    fn new(id: &str, seq: PyRef<'_, Protein>, desc: Option<&str>) -> Self {
        let mut record = SeqRecord::new(id, seq.inner.clone());
        if let Some(desc) = desc {
            record = record.with_desc(desc);
        }
        Self { inner: record }
    }

    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }

    #[getter]
    fn description(&self) -> Option<&str> {
        self.inner.desc.as_deref()
    }

    #[getter]
    fn seq(&self) -> Protein {
        Protein {
            inner: self.inner.seq.clone(),
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let seq_repr = seq_shared::seq_repr(self.inner.seq.as_bytes(), "Protein");
        let desc = match self.inner.desc.as_deref() {
            Some(desc) => format!("{desc:?}"),
            None => "None".to_string(),
        };
        Ok(format!(
            "ProteinRecord(id={:?}, description={desc}, seq={seq_repr})",
            self.inner.id
        ))
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ProteinRecord>()?;
    Ok(())
}
