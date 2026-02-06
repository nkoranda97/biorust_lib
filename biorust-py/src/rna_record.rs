use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule};

use crate::feature;
use crate::feature::SeqFeature;
use crate::rna::RNA;
use crate::seq_shared;
use biorust_core::seq::record::SeqRecord;
use biorust_core::seq::rna::RnaSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
pub struct RNARecord {
    pub(crate) inner: SeqRecord<RnaSeq>,
}

#[pymethods]
impl RNARecord {
    #[new]
    #[pyo3(signature = (id, seq, desc=None, features=None, annotations=None))]
    fn new(
        id: &str,
        seq: PyRef<'_, RNA>,
        desc: Option<&str>,
        features: Option<&Bound<'_, PyAny>>,
        annotations: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut record = SeqRecord::new(id, seq.inner.clone());
        if let Some(desc) = desc {
            record = record.with_desc(desc);
        }
        if let Some(features) = features {
            record = record.with_features(feature::extract_features(features)?);
        }
        if let Some(annotations) = annotations {
            record =
                record.with_annotations(feature::extract_str_list_map(annotations, "annotations")?);
        }
        Ok(Self { inner: record })
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
    fn seq(&self) -> RNA {
        RNA {
            inner: self.inner.seq.clone(),
        }
    }

    #[getter]
    fn features(&self, py: Python<'_>) -> PyResult<Vec<Py<SeqFeature>>> {
        feature::features_to_pylist(py, &self.inner.features)
    }

    #[getter]
    fn annotations(&self, py: Python<'_>) -> PyResult<PyObject> {
        feature::map_to_pydict(py, &self.inner.annotations)
    }

    fn __repr__(&self) -> PyResult<String> {
        let seq_repr = seq_shared::seq_repr(self.inner.seq.as_bytes(), "RNA");
        let desc = match self.inner.desc.as_deref() {
            Some(desc) => format!("{desc:?}"),
            None => "None".to_string(),
        };
        Ok(format!(
            "RNARecord(id={:?}, description={desc}, seq={seq_repr})",
            self.inner.id
        ))
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RNARecord>()?;
    Ok(())
}
