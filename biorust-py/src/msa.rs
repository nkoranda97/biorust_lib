use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyModule, PySlice, PyTuple};

use crate::gapped_dna::GappedDNA;
use crate::gapped_protein::GappedProtein;
use biorust_core::seq::gapped_dna::GappedDnaSeq;
use biorust_core::seq::gapped_protein::GappedProteinSeq;

#[inline]
fn is_gap(b: u8) -> bool {
    b == b'-' || b == b'.'
}

#[allow(clippy::upper_case_acronyms)]
#[pyclass(frozen)]
pub struct AlignmentDNA {
    ids: Vec<Box<str>>,
    seqs: Vec<GappedDnaSeq>,
    width: usize,
}

impl AlignmentDNA {
    pub(crate) fn seqs_ref(&self) -> &[GappedDnaSeq] {
        &self.seqs
    }

    pub(crate) fn labels_cloned(&self) -> Vec<Box<str>> {
        self.ids.clone()
    }
}

#[pymethods]
impl AlignmentDNA {
    #[new]
    fn new(records: &Bound<'_, PyList>) -> PyResult<Self> {
        let mut ids = Vec::with_capacity(records.len());
        let mut seqs = Vec::with_capacity(records.len());

        for item in records.iter() {
            let tuple = item.downcast::<PyTuple>().map_err(|_| {
                PyValueError::new_err("AlignmentDNA expects a list of (str, GappedDNA) tuples")
            })?;

            if tuple.len() != 2 {
                return Err(PyValueError::new_err(
                    "each tuple must have exactly 2 elements: (id, seq)",
                ));
            }

            let id: String = tuple.get_item(0)?.extract()?;
            let gapped: PyRef<'_, GappedDNA> = tuple.get_item(1)?.extract()?;

            ids.push(id.into_boxed_str());
            seqs.push(gapped.inner.clone());
        }

        if seqs.is_empty() {
            return Err(PyValueError::new_err(
                "AlignmentDNA requires at least one sequence",
            ));
        }

        let width = seqs[0].len();
        for (i, s) in seqs.iter().enumerate() {
            if s.len() != width {
                return Err(PyValueError::new_err(format!(
                    "sequence {} has length {} but expected {} (all sequences must be equal length)",
                    i,
                    s.len(),
                    width,
                )));
            }
        }

        Ok(Self { ids, seqs, width })
    }

    #[getter]
    fn width(&self) -> usize {
        self.width
    }

    fn __len__(&self) -> usize {
        self.ids.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if let Ok(slice) = index.downcast::<PySlice>() {
            let n = self.ids.len() as isize;
            let idx = slice.indices(n)?;
            let (start, stop, step) = (idx.start, idx.stop, idx.step);

            let mut new_ids = Vec::new();
            let mut new_seqs = Vec::new();

            if step > 0 {
                let mut i = start;
                while i < stop {
                    new_ids.push(self.ids[i as usize].clone());
                    new_seqs.push(self.seqs[i as usize].clone());
                    i += step;
                }
            } else {
                let mut i = start;
                while i > stop {
                    new_ids.push(self.ids[i as usize].clone());
                    new_seqs.push(self.seqs[i as usize].clone());
                    i += step;
                }
            }

            let width = new_seqs.first().map(|s| s.len()).unwrap_or(0);
            let obj = AlignmentDNA {
                ids: new_ids,
                seqs: new_seqs,
                width,
            };
            return Ok(obj.into_py(py));
        }

        let i: isize = index
            .extract()
            .map_err(|_| PyTypeError::new_err("index must be int or slice"))?;

        let n = self.ids.len() as isize;
        let idx = if i < 0 { i + n } else { i };
        if idx < 0 || idx >= n {
            return Err(PyIndexError::new_err("index out of range"));
        }
        let idx = idx as usize;

        let id = self.ids[idx].as_ref();
        let gapped = GappedDNA {
            inner: self.seqs[idx].clone(),
        };

        let tuple = PyTuple::new_bound(py, &[id.to_object(py), gapped.into_py(py)]);
        Ok(tuple.to_object(py))
    }

    fn ids(&self) -> Vec<String> {
        self.ids.iter().map(|s| s.to_string()).collect()
    }

    fn seqs(&self) -> Vec<GappedDNA> {
        self.seqs
            .iter()
            .map(|s| GappedDNA { inner: s.clone() })
            .collect()
    }

    /// Return all aligned sequences as a list of strings.
    fn aligned_strings(&self) -> Vec<String> {
        self.seqs
            .iter()
            .map(|s| {
                std::str::from_utf8(s.as_bytes())
                    .unwrap_or("<bytes>")
                    .to_string()
            })
            .collect()
    }

    /// Build a Clustal-style alignment diagram.
    ///
    /// Each sequence is labelled with its ID (left-padded), followed by a
    /// conservation line where ``*`` marks fully-conserved columns.
    fn alignment_diagram(&self) -> String {
        if self.seqs.is_empty() {
            return String::new();
        }
        let pad = self.ids.iter().map(|id| id.len()).max().unwrap_or(0);
        let mut lines = Vec::with_capacity(self.ids.len() + 1);

        for (id, seq) in self.ids.iter().zip(&self.seqs) {
            let s = std::str::from_utf8(seq.as_bytes()).unwrap_or("<bytes>");
            lines.push(format!("{:>pad$}  {}", id, s, pad = pad));
        }

        // Conservation line: '*' if all bases in column are identical
        // (case-insensitive), ' ' otherwise.
        let mut conservation = String::with_capacity(self.width);
        for col in 0..self.width {
            let first = self.seqs[0].as_bytes()[col].to_ascii_uppercase();
            let conserved = self.seqs[1..]
                .iter()
                .all(|s| s.as_bytes()[col].to_ascii_uppercase() == first && !is_gap(first));
            conservation.push(if conserved { '*' } else { ' ' });
        }
        lines.push(format!("{:>pad$}  {}", "", conservation, pad = pad));

        lines.join("\n")
    }

    fn __repr__(&self) -> String {
        format!("AlignmentDNA(n={}, width={})", self.ids.len(), self.width)
    }

    fn __str__(&self) -> String {
        self.alignment_diagram()
    }
}

#[pyclass(frozen)]
pub struct AlignmentProtein {
    ids: Vec<Box<str>>,
    seqs: Vec<GappedProteinSeq>,
    width: usize,
}

impl AlignmentProtein {
    pub(crate) fn seqs_ref(&self) -> &[GappedProteinSeq] {
        &self.seqs
    }

    pub(crate) fn labels_cloned(&self) -> Vec<Box<str>> {
        self.ids.clone()
    }
}

#[pymethods]
impl AlignmentProtein {
    #[new]
    fn new(records: &Bound<'_, PyList>) -> PyResult<Self> {
        let mut ids = Vec::with_capacity(records.len());
        let mut seqs = Vec::with_capacity(records.len());

        for item in records.iter() {
            let tuple = item.downcast::<PyTuple>().map_err(|_| {
                PyValueError::new_err(
                    "AlignmentProtein expects a list of (str, GappedProtein) tuples",
                )
            })?;

            if tuple.len() != 2 {
                return Err(PyValueError::new_err(
                    "each tuple must have exactly 2 elements: (id, seq)",
                ));
            }

            let id: String = tuple.get_item(0)?.extract()?;
            let gapped: PyRef<'_, GappedProtein> = tuple.get_item(1)?.extract()?;

            ids.push(id.into_boxed_str());
            seqs.push(gapped.inner.clone());
        }

        if seqs.is_empty() {
            return Err(PyValueError::new_err(
                "AlignmentProtein requires at least one sequence",
            ));
        }

        let width = seqs[0].len();
        for (i, s) in seqs.iter().enumerate() {
            if s.len() != width {
                return Err(PyValueError::new_err(format!(
                    "sequence {} has length {} but expected {} (all sequences must be equal length)",
                    i,
                    s.len(),
                    width,
                )));
            }
        }

        Ok(Self { ids, seqs, width })
    }

    #[getter]
    fn width(&self) -> usize {
        self.width
    }

    fn __len__(&self) -> usize {
        self.ids.len()
    }

    fn __getitem__(&self, py: Python<'_>, index: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if let Ok(slice) = index.downcast::<PySlice>() {
            let n = self.ids.len() as isize;
            let idx = slice.indices(n)?;
            let (start, stop, step) = (idx.start, idx.stop, idx.step);

            let mut new_ids = Vec::new();
            let mut new_seqs = Vec::new();

            if step > 0 {
                let mut i = start;
                while i < stop {
                    new_ids.push(self.ids[i as usize].clone());
                    new_seqs.push(self.seqs[i as usize].clone());
                    i += step;
                }
            } else {
                let mut i = start;
                while i > stop {
                    new_ids.push(self.ids[i as usize].clone());
                    new_seqs.push(self.seqs[i as usize].clone());
                    i += step;
                }
            }

            let width = new_seqs.first().map(|s| s.len()).unwrap_or(0);
            let obj = AlignmentProtein {
                ids: new_ids,
                seqs: new_seqs,
                width,
            };
            return Ok(obj.into_py(py));
        }

        let i: isize = index
            .extract()
            .map_err(|_| PyTypeError::new_err("index must be int or slice"))?;

        let n = self.ids.len() as isize;
        let idx = if i < 0 { i + n } else { i };
        if idx < 0 || idx >= n {
            return Err(PyIndexError::new_err("index out of range"));
        }
        let idx = idx as usize;

        let id = self.ids[idx].as_ref();
        let gapped = GappedProtein {
            inner: self.seqs[idx].clone(),
        };

        let tuple = PyTuple::new_bound(py, &[id.to_object(py), gapped.into_py(py)]);
        Ok(tuple.to_object(py))
    }

    fn ids(&self) -> Vec<String> {
        self.ids.iter().map(|s| s.to_string()).collect()
    }

    fn seqs(&self) -> Vec<GappedProtein> {
        self.seqs
            .iter()
            .map(|s| GappedProtein { inner: s.clone() })
            .collect()
    }

    fn aligned_strings(&self) -> Vec<String> {
        self.seqs
            .iter()
            .map(|s| {
                std::str::from_utf8(s.as_bytes())
                    .unwrap_or("<bytes>")
                    .to_string()
            })
            .collect()
    }

    fn alignment_diagram(&self) -> String {
        if self.seqs.is_empty() {
            return String::new();
        }
        let pad = self.ids.iter().map(|id| id.len()).max().unwrap_or(0);
        let mut lines = Vec::with_capacity(self.ids.len() + 1);

        for (id, seq) in self.ids.iter().zip(&self.seqs) {
            let s = std::str::from_utf8(seq.as_bytes()).unwrap_or("<bytes>");
            lines.push(format!("{:>pad$}  {}", id, s, pad = pad));
        }

        let mut conservation = String::with_capacity(self.width);
        for col in 0..self.width {
            let first = self.seqs[0].as_bytes()[col].to_ascii_uppercase();
            let conserved = self.seqs[1..]
                .iter()
                .all(|s| s.as_bytes()[col].to_ascii_uppercase() == first && !is_gap(first));
            conservation.push(if conserved { '*' } else { ' ' });
        }
        lines.push(format!("{:>pad$}  {}", "", conservation, pad = pad));

        lines.join("\n")
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentProtein(n={}, width={})",
            self.ids.len(),
            self.width
        )
    }

    fn __str__(&self) -> String {
        self.alignment_diagram()
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<AlignmentDNA>()?;
    m.add_class::<AlignmentProtein>()?;
    Ok(())
}
