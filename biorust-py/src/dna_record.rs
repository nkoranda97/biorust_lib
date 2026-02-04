use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::dna::DNA;
use crate::seq_shared;
use biorust_core::seq::dna::DnaSeq;
use biorust_core::seq::record::SeqRecord;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
pub struct DNARecord {
    pub(crate) inner: SeqRecord<DnaSeq>,
}

#[pymethods]
impl DNARecord {
    #[new]
    #[pyo3(signature = (id, seq, desc=None))]
    fn new(id: &str, seq: PyRef<'_, DNA>, desc: Option<&str>) -> Self {
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
    fn seq(&self) -> DNA {
        DNA {
            inner: self.inner.seq.clone(),
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let seq_repr = seq_shared::seq_repr(self.inner.seq.as_bytes(), "DNA");
        let desc = match self.inner.desc.as_deref() {
            Some(desc) => format!("{desc:?}"),
            None => "None".to_string(),
        };
        Ok(format!(
            "DNARecord(id={:?}, description={desc}, seq={seq_repr})",
            self.inner.id
        ))
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DNARecord>()?;
    Ok(())
}
