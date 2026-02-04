use pyo3::prelude::*;
use pyo3::types::PyModule;

use biorust_core::io::SkippedRecord as CoreSkippedRecord;

#[pyclass(frozen)]
#[derive(Clone)]
pub struct SkippedRecord {
    pub(crate) row: usize,
    pub(crate) id: Option<String>,
    pub(crate) column: String,
    pub(crate) message: String,
}

#[pymethods]
impl SkippedRecord {
    #[getter]
    fn row(&self) -> usize {
        self.row
    }

    #[getter]
    fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }

    #[getter]
    fn column(&self) -> &str {
        &self.column
    }

    #[getter]
    fn message(&self) -> &str {
        &self.message
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "SkippedRecord(row={}, id={:?}, column={:?}, message={:?})",
            self.row, self.id, self.column, self.message
        ))
    }
}

impl From<CoreSkippedRecord> for SkippedRecord {
    fn from(value: CoreSkippedRecord) -> Self {
        Self {
            row: value.row,
            id: value.id.map(|s| s.to_string()),
            column: value.column.to_string(),
            message: value.message.to_string(),
        }
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SkippedRecord>()?;
    Ok(())
}
