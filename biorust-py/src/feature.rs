use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyList, PyModule};
use std::collections::HashMap;

use biorust_core::seq::feature::{
    FeatureLocation as CoreFeatureLocation, Qualifiers, SeqFeature as CoreSeqFeature,
};

#[pyclass(frozen)]
pub struct FeatureLocation {
    pub(crate) inner: CoreFeatureLocation,
}

#[pymethods]
impl FeatureLocation {
    #[new]
    #[pyo3(signature = (start, end, strand=None))]
    fn new(start: usize, end: usize, strand: Option<i8>) -> PyResult<Self> {
        let inner = CoreFeatureLocation::new(start, end, strand)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    #[getter]
    fn start(&self) -> usize {
        self.inner.start()
    }

    #[getter]
    fn end(&self) -> usize {
        self.inner.end()
    }

    #[getter]
    fn strand(&self) -> Option<i8> {
        self.inner.strand()
    }

    fn __repr__(&self) -> String {
        format!(
            "FeatureLocation(start={}, end={}, strand={:?})",
            self.inner.start(),
            self.inner.end(),
            self.inner.strand()
        )
    }
}

#[pyclass(frozen)]
pub struct SeqFeature {
    pub(crate) inner: CoreSeqFeature,
}

#[pymethods]
impl SeqFeature {
    #[new]
    #[pyo3(signature = (feature_type, location, qualifiers=None))]
    fn new(
        feature_type: &str,
        location: PyRef<'_, FeatureLocation>,
        qualifiers: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut inner = CoreSeqFeature::new(feature_type, location.inner.clone())
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        if let Some(qualifiers) = qualifiers {
            let parsed = extract_str_list_map(qualifiers, "qualifiers")?;
            inner = inner.with_qualifiers(parsed);
        }
        Ok(Self { inner })
    }

    #[getter]
    fn feature_type(&self) -> &str {
        self.inner.feature_type()
    }

    #[getter]
    fn location(&self) -> FeatureLocation {
        FeatureLocation {
            inner: self.inner.location().clone(),
        }
    }

    #[getter]
    fn qualifiers(&self, py: Python<'_>) -> PyResult<PyObject> {
        map_to_pydict(py, self.inner.qualifiers())
    }

    fn __repr__(&self) -> String {
        format!(
            "SeqFeature(type={:?}, location={})",
            self.inner.feature_type(),
            FeatureLocation {
                inner: self.inner.location().clone()
            }
            .__repr__()
        )
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FeatureLocation>()?;
    m.add_class::<SeqFeature>()?;
    Ok(())
}

pub(crate) fn extract_str_list_map(
    obj: &Bound<'_, PyAny>,
    name: &str,
) -> PyResult<HashMap<Box<str>, Vec<Box<str>>>> {
    let dict = obj
        .downcast::<PyDict>()
        .map_err(|_| PyTypeError::new_err(format!("{name} must be a dict[str, list[str]]")))?;
    let mut out = HashMap::new();
    for (key, value) in dict.iter() {
        let key: String = key
            .extract()
            .map_err(|_| PyTypeError::new_err(format!("{name} keys must be str")))?;
        let values: Vec<String> = value
            .extract()
            .map_err(|_| PyTypeError::new_err(format!("{name} values must be list[str]")))?;
        let boxed_values = values
            .into_iter()
            .map(|v| v.into_boxed_str())
            .collect::<Vec<_>>();
        out.insert(key.into_boxed_str(), boxed_values);
    }
    Ok(out)
}

pub(crate) fn map_to_pydict(py: Python<'_>, map: &Qualifiers) -> PyResult<PyObject> {
    let dict = PyDict::new_bound(py);
    for (key, values) in map {
        let list = PyList::new_bound(py, values.iter().map(|v| v.as_ref()));
        dict.set_item(key.as_ref(), list)?;
    }
    Ok(dict.to_object(py))
}

pub(crate) fn extract_features(obj: &Bound<'_, PyAny>) -> PyResult<Vec<CoreSeqFeature>> {
    let mut out = Vec::new();
    for item in obj.iter()? {
        let item = item?;
        let feature = item
            .extract::<PyRef<'_, SeqFeature>>()
            .map_err(|_| PyTypeError::new_err("features must be SeqFeature objects"))?;
        out.push(feature.inner.clone());
    }
    Ok(out)
}

pub(crate) fn features_to_pylist(
    py: Python<'_>,
    features: &[CoreSeqFeature],
) -> PyResult<Vec<Py<SeqFeature>>> {
    let mut out = Vec::with_capacity(features.len());
    for feature in features {
        out.push(Py::new(
            py,
            SeqFeature {
                inner: feature.clone(),
            },
        )?);
    }
    Ok(out)
}
