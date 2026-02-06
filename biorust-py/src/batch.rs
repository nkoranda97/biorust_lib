#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList, PySlice};

use crate::dna::DNA;
use crate::protein::Protein;
use crate::rna::RNA;
use biorust_core::seq::batch::SeqBatch;
use biorust_core::seq::dna::DnaSeq;
use biorust_core::seq::protein::ProteinSeq;
use biorust_core::seq::rna::RnaSeq;

#[allow(clippy::upper_case_acronyms)]
#[pyclass]
pub struct DNABatch {
    pub(crate) inner: SeqBatch<DnaSeq>,
}

#[allow(clippy::upper_case_acronyms)]
#[pyclass]
pub struct ProteinBatch {
    pub(crate) inner: SeqBatch<ProteinSeq>,
}

#[allow(clippy::upper_case_acronyms)]
#[pyclass]
pub struct RNABatch {
    pub(crate) inner: SeqBatch<RnaSeq>,
}

fn collect_dna_seqs(obj: &Bound<'_, PyAny>) -> PyResult<Vec<DnaSeq>> {
    if let Ok(batch) = obj.extract::<PyRef<'_, DNABatch>>() {
        return Ok(batch.inner.as_slice().to_vec());
    }

    let mut out = Vec::new();
    for item in obj.iter()? {
        let item = item?;
        let dna = item
            .extract::<PyRef<'_, DNA>>()
            .map_err(|_| PyTypeError::new_err("DNABatch expects DNA objects only"))?;
        out.push(dna.inner.clone());
    }
    Ok(out)
}

fn collect_rna_seqs(obj: &Bound<'_, PyAny>) -> PyResult<Vec<RnaSeq>> {
    if let Ok(batch) = obj.extract::<PyRef<'_, RNABatch>>() {
        return Ok(batch.inner.as_slice().to_vec());
    }

    let mut out = Vec::new();
    for item in obj.iter()? {
        let item = item?;
        let rna = item
            .extract::<PyRef<'_, RNA>>()
            .map_err(|_| PyTypeError::new_err("RNABatch expects RNA objects only"))?;
        out.push(rna.inner.clone());
    }
    Ok(out)
}

fn collect_protein_seqs(obj: &Bound<'_, PyAny>) -> PyResult<Vec<ProteinSeq>> {
    if let Ok(batch) = obj.extract::<PyRef<'_, ProteinBatch>>() {
        return Ok(batch.inner.as_slice().to_vec());
    }

    let mut out = Vec::new();
    for item in obj.iter()? {
        let item = item?;
        let protein = item
            .extract::<PyRef<'_, Protein>>()
            .map_err(|_| PyTypeError::new_err("ProteinBatch expects Protein objects only"))?;
        out.push(protein.inner.clone());
    }
    Ok(out)
}

fn normalize_slice(
    len: usize,
    start: Option<isize>,
    stop: Option<isize>,
    step: isize,
) -> PyResult<(usize, usize, usize)> {
    if step <= 0 {
        return Err(PyValueError::new_err("step must be >= 1"));
    }

    let n = len as isize;
    let mut s = start.unwrap_or(0);
    let mut e = stop.unwrap_or(n);

    if s < 0 {
        s += n;
    }
    if e < 0 {
        e += n;
    }

    s = s.clamp(0, n);
    e = e.clamp(0, n);

    Ok((s as usize, e as usize, step as usize))
}

fn collect_take_indices(obj: &Bound<'_, PyAny>, len: usize) -> PyResult<Vec<usize>> {
    let iter = obj
        .iter()
        .map_err(|_| PyTypeError::new_err("idxs must be an iterable of ints"))?;
    let mut out = Vec::new();
    let n = len as isize;

    for item in iter {
        let idx: isize = item
            .map_err(|_| PyTypeError::new_err("idxs must be an iterable of ints"))?
            .extract()
            .map_err(|_| PyTypeError::new_err("idxs must be an iterable of ints"))?;
        let idx = if idx < 0 { idx + n } else { idx };
        if idx < 0 || idx >= n {
            return Err(PyIndexError::new_err("index out of range"));
        }
        out.push(idx as usize);
    }

    Ok(out)
}

#[pymethods]
impl DNABatch {
    #[new]
    fn new(seqs: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: SeqBatch::new(collect_dna_seqs(seqs)?),
        })
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        if let Ok(slice) = index.downcast::<PySlice>() {
            let idx = slice.indices(self.inner.len() as isize)?;
            let (start, stop, step) = (idx.start, idx.stop, idx.step);
            let mut out = Vec::new();

            if step > 0 {
                let mut i = start;
                while i < stop {
                    out.push(self.inner[i as usize].clone());
                    i += step;
                }
            } else {
                let mut i = start;
                while i > stop {
                    out.push(self.inner[i as usize].clone());
                    i += step;
                }
            }

            let batch = DNABatch {
                inner: SeqBatch::new(out),
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

        Ok(Py::new(
            py,
            DNA {
                inner: self.inner[i as usize].clone(),
            },
        )?
        .to_object(py))
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let list = self.to_list(py)?;
        list.call_method0("__iter__")
    }

    fn to_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let mut items = Vec::with_capacity(self.inner.len());
        for seq in self.inner.as_slice() {
            items.push(Py::new(py, DNA { inner: seq.clone() })?);
        }
        Ok(PyList::new_bound(py, items))
    }

    fn lengths(&self) -> Vec<usize> {
        self.inner.lengths()
    }

    fn copy(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }

    #[pyo3(signature = (start=None, stop=None, step=1))]
    fn slice(&self, start: Option<isize>, stop: Option<isize>, step: isize) -> PyResult<Self> {
        let (start, stop, step) = normalize_slice(self.inner.len(), start, stop, step)?;
        Ok(Self {
            inner: self.inner.slice(start, stop, step),
        })
    }

    fn take(&self, idxs: &Bound<'_, PyAny>) -> PyResult<Self> {
        let idxs = collect_take_indices(idxs, self.inner.len())?;
        let out = self
            .inner
            .take(&idxs)
            .map_err(|err| PyIndexError::new_err(err.to_string()))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (min_len=None, max_len=None, inplace=false))]
    fn filter_by_len(
        &mut self,
        py: Python<'_>,
        min_len: Option<usize>,
        max_len: Option<usize>,
        inplace: bool,
    ) -> PyResult<PyObject> {
        let filtered = self.inner.filter_by_len(min_len, max_len);
        if inplace {
            self.inner = filtered;
            return Ok(py.None());
        }
        let out = DNABatch { inner: filtered };
        Ok(Py::new(py, out)?.to_object(py))
    }

    fn concat(&self, py: Python<'_>) -> PyResult<PyObject> {
        let seq = self
            .inner
            .concat_all()
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Py::new(py, DNA { inner: seq })?.to_object(py))
    }

    fn count(&self, needle: &Bound<'_, PyAny>) -> PyResult<Vec<usize>> {
        let dna = needle
            .extract::<PyRef<'_, DNA>>()
            .map_err(|_| PyTypeError::new_err("DNABatch.count expects a DNA object"))?;
        self.inner
            .count(&dna.inner)
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }

    fn contains(&self, needle: &Bound<'_, PyAny>) -> PyResult<Vec<bool>> {
        let dna = needle
            .extract::<PyRef<'_, DNA>>()
            .map_err(|_| PyTypeError::new_err("DNABatch.contains expects a DNA object"))?;
        self.inner
            .contains(&dna.inner)
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }

    fn append(&mut self, seq: &Bound<'_, PyAny>) -> PyResult<()> {
        let dna = seq
            .extract::<PyRef<'_, DNA>>()
            .map_err(|_| PyTypeError::new_err("DNABatch expects DNA objects only"))?;
        self.inner.push(dna.inner.clone());
        Ok(())
    }

    fn extend(&mut self, seqs: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = collect_dna_seqs(seqs)?;
        self.inner.extend(out);
        Ok(())
    }

    fn clear(&mut self) {
        self.inner.clear();
    }

    fn reserve(&mut self, additional: usize) {
        self.inner.reserve(additional);
    }

    fn pop(&mut self, py: Python<'_>) -> PyResult<PyObject> {
        match self.inner.pop() {
            Some(seq) => Ok(Py::new(py, DNA { inner: seq })?.to_object(py)),
            None => Err(PyIndexError::new_err("pop from empty batch")),
        }
    }

    fn truncate(&mut self, len: usize) {
        self.inner.truncate(len);
    }

    fn __iadd__(mut slf: PyRefMut<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = collect_dna_seqs(other)?;
        slf.inner.extend(out);
        Ok(())
    }

    fn __imul__(mut slf: PyRefMut<'_, Self>, n: isize) -> PyResult<()> {
        if n <= 0 {
            slf.inner.clear();
            return Ok(());
        }
        if n == 1 {
            return Ok(());
        }
        let n = n as usize;
        let orig: Vec<DnaSeq> = slf.inner.as_slice().to_vec();
        slf.inner.clear();
        slf.inner.reserve(orig.len() * n);
        for _ in 0..n {
            slf.inner.extend(orig.iter().cloned());
        }
        Ok(())
    }

    #[pyo3(signature = (inplace=false))]
    fn reverse_complements(&mut self, py: Python<'_>, inplace: bool) -> PyResult<PyObject> {
        if inplace {
            self.inner.reverse_complements_in_place();
            return Ok(py.None());
        }

        let out = DNABatch {
            inner: self.inner.reverse_complements(),
        };
        Ok(Py::new(py, out)?.to_object(py))
    }
}

#[pymethods]
impl RNABatch {
    #[new]
    fn new(seqs: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: SeqBatch::new(collect_rna_seqs(seqs)?),
        })
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        if let Ok(slice) = index.downcast::<PySlice>() {
            let idx = slice.indices(self.inner.len() as isize)?;
            let (start, stop, step) = (idx.start, idx.stop, idx.step);
            let mut out = Vec::new();

            if step > 0 {
                let mut i = start;
                while i < stop {
                    out.push(self.inner[i as usize].clone());
                    i += step;
                }
            } else {
                let mut i = start;
                while i > stop {
                    out.push(self.inner[i as usize].clone());
                    i += step;
                }
            }

            let batch = RNABatch {
                inner: SeqBatch::new(out),
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

        Ok(Py::new(
            py,
            RNA {
                inner: self.inner[i as usize].clone(),
            },
        )?
        .to_object(py))
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let list = self.to_list(py)?;
        list.call_method0("__iter__")
    }

    fn to_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let mut items = Vec::with_capacity(self.inner.len());
        for seq in self.inner.as_slice() {
            items.push(Py::new(py, RNA { inner: seq.clone() })?);
        }
        Ok(PyList::new_bound(py, items))
    }

    fn lengths(&self) -> Vec<usize> {
        self.inner.lengths()
    }

    fn copy(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }

    #[pyo3(signature = (start=None, stop=None, step=1))]
    fn slice(&self, start: Option<isize>, stop: Option<isize>, step: isize) -> PyResult<Self> {
        let (start, stop, step) = normalize_slice(self.inner.len(), start, stop, step)?;
        Ok(Self {
            inner: self.inner.slice(start, stop, step),
        })
    }

    fn take(&self, idxs: &Bound<'_, PyAny>) -> PyResult<Self> {
        let idxs = collect_take_indices(idxs, self.inner.len())?;
        let out = self
            .inner
            .take(&idxs)
            .map_err(|err| PyIndexError::new_err(err.to_string()))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (min_len=None, max_len=None, inplace=false))]
    fn filter_by_len(
        &mut self,
        py: Python<'_>,
        min_len: Option<usize>,
        max_len: Option<usize>,
        inplace: bool,
    ) -> PyResult<PyObject> {
        let filtered = self.inner.filter_by_len(min_len, max_len);
        if inplace {
            self.inner = filtered;
            return Ok(py.None());
        }
        let out = RNABatch { inner: filtered };
        Ok(Py::new(py, out)?.to_object(py))
    }

    fn concat(&self, py: Python<'_>) -> PyResult<PyObject> {
        let seq = self
            .inner
            .concat_all()
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Py::new(py, RNA { inner: seq })?.to_object(py))
    }

    fn append(&mut self, seq: &Bound<'_, PyAny>) -> PyResult<()> {
        let rna = seq
            .extract::<PyRef<'_, RNA>>()
            .map_err(|_| PyTypeError::new_err("RNABatch expects RNA objects only"))?;
        self.inner.push(rna.inner.clone());
        Ok(())
    }

    fn extend(&mut self, seqs: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = collect_rna_seqs(seqs)?;
        self.inner.extend(out);
        Ok(())
    }

    fn clear(&mut self) {
        self.inner.clear();
    }

    fn reserve(&mut self, additional: usize) {
        self.inner.reserve(additional);
    }

    fn pop(&mut self, py: Python<'_>) -> PyResult<PyObject> {
        match self.inner.pop() {
            Some(seq) => Ok(Py::new(py, RNA { inner: seq })?.to_object(py)),
            None => Err(PyIndexError::new_err("pop from empty batch")),
        }
    }

    fn truncate(&mut self, len: usize) {
        self.inner.truncate(len);
    }

    fn __iadd__(mut slf: PyRefMut<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = collect_rna_seqs(other)?;
        slf.inner.extend(out);
        Ok(())
    }

    fn __imul__(mut slf: PyRefMut<'_, Self>, n: isize) -> PyResult<()> {
        if n <= 0 {
            slf.inner.clear();
            return Ok(());
        }
        if n == 1 {
            return Ok(());
        }
        let n = n as usize;
        let orig: Vec<RnaSeq> = slf.inner.as_slice().to_vec();
        slf.inner.clear();
        slf.inner.reserve(orig.len() * n);
        for _ in 0..n {
            slf.inner.extend(orig.iter().cloned());
        }
        Ok(())
    }

    #[pyo3(signature = (inplace=false))]
    fn reverse_complements(&mut self, py: Python<'_>, inplace: bool) -> PyResult<PyObject> {
        if inplace {
            self.inner.reverse_complements_in_place();
            return Ok(py.None());
        }

        let out = RNABatch {
            inner: self.inner.reverse_complements(),
        };
        Ok(Py::new(py, out)?.to_object(py))
    }
}

#[pymethods]
impl ProteinBatch {
    #[new]
    fn new(seqs: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: SeqBatch::new(collect_protein_seqs(seqs)?),
        })
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __getitem__<'py>(&self, py: Python<'py>, index: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        if let Ok(slice) = index.downcast::<PySlice>() {
            let idx = slice.indices(self.inner.len() as isize)?;
            let (start, stop, step) = (idx.start, idx.stop, idx.step);
            let mut out = Vec::new();

            if step > 0 {
                let mut i = start;
                while i < stop {
                    out.push(self.inner[i as usize].clone());
                    i += step;
                }
            } else {
                let mut i = start;
                while i > stop {
                    out.push(self.inner[i as usize].clone());
                    i += step;
                }
            }

            let batch = ProteinBatch {
                inner: SeqBatch::new(out),
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

        Ok(Py::new(
            py,
            Protein {
                inner: self.inner[i as usize].clone(),
            },
        )?
        .to_object(py))
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let list = self.to_list(py)?;
        list.call_method0("__iter__")
    }

    fn to_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let mut items = Vec::with_capacity(self.inner.len());
        for seq in self.inner.as_slice() {
            items.push(Py::new(py, Protein { inner: seq.clone() })?);
        }
        Ok(PyList::new_bound(py, items))
    }

    fn lengths(&self) -> Vec<usize> {
        self.inner.lengths()
    }

    fn copy(&self) -> Self {
        Self {
            inner: self.inner.clone(),
        }
    }

    #[pyo3(signature = (start=None, stop=None, step=1))]
    fn slice(&self, start: Option<isize>, stop: Option<isize>, step: isize) -> PyResult<Self> {
        let (start, stop, step) = normalize_slice(self.inner.len(), start, stop, step)?;
        Ok(Self {
            inner: self.inner.slice(start, stop, step),
        })
    }

    fn take(&self, idxs: &Bound<'_, PyAny>) -> PyResult<Self> {
        let idxs = collect_take_indices(idxs, self.inner.len())?;
        let out = self
            .inner
            .take(&idxs)
            .map_err(|err| PyIndexError::new_err(err.to_string()))?;
        Ok(Self { inner: out })
    }

    #[pyo3(signature = (min_len=None, max_len=None, inplace=false))]
    fn filter_by_len(
        &mut self,
        py: Python<'_>,
        min_len: Option<usize>,
        max_len: Option<usize>,
        inplace: bool,
    ) -> PyResult<PyObject> {
        let filtered = self.inner.filter_by_len(min_len, max_len);
        if inplace {
            self.inner = filtered;
            return Ok(py.None());
        }
        let out = ProteinBatch { inner: filtered };
        Ok(Py::new(py, out)?.to_object(py))
    }

    fn concat(&self, py: Python<'_>) -> PyResult<PyObject> {
        let seq = self
            .inner
            .concat_all()
            .map_err(|err| PyValueError::new_err(err.to_string()))?;
        Ok(Py::new(py, Protein { inner: seq })?.to_object(py))
    }

    fn count(&self, needle: &Bound<'_, PyAny>) -> PyResult<Vec<usize>> {
        let protein = needle
            .extract::<PyRef<'_, Protein>>()
            .map_err(|_| PyTypeError::new_err("ProteinBatch.count expects a Protein object"))?;
        self.inner
            .count(&protein.inner)
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }

    fn contains(&self, needle: &Bound<'_, PyAny>) -> PyResult<Vec<bool>> {
        let protein = needle
            .extract::<PyRef<'_, Protein>>()
            .map_err(|_| PyTypeError::new_err("ProteinBatch.contains expects a Protein object"))?;
        self.inner
            .contains(&protein.inner)
            .map_err(|err| PyValueError::new_err(err.to_string()))
    }

    fn append(&mut self, seq: &Bound<'_, PyAny>) -> PyResult<()> {
        let protein = seq
            .extract::<PyRef<'_, Protein>>()
            .map_err(|_| PyTypeError::new_err("ProteinBatch expects Protein objects only"))?;
        self.inner.push(protein.inner.clone());
        Ok(())
    }

    fn extend(&mut self, seqs: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = collect_protein_seqs(seqs)?;
        self.inner.extend(out);
        Ok(())
    }

    fn clear(&mut self) {
        self.inner.clear();
    }

    fn reserve(&mut self, additional: usize) {
        self.inner.reserve(additional);
    }

    fn pop(&mut self, py: Python<'_>) -> PyResult<PyObject> {
        match self.inner.pop() {
            Some(seq) => Ok(Py::new(py, Protein { inner: seq })?.to_object(py)),
            None => Err(PyIndexError::new_err("pop from empty batch")),
        }
    }

    fn truncate(&mut self, len: usize) {
        self.inner.truncate(len);
    }

    fn __iadd__(mut slf: PyRefMut<'_, Self>, other: &Bound<'_, PyAny>) -> PyResult<()> {
        let out = collect_protein_seqs(other)?;
        slf.inner.extend(out);
        Ok(())
    }

    fn __imul__(mut slf: PyRefMut<'_, Self>, n: isize) -> PyResult<()> {
        if n <= 0 {
            slf.inner.clear();
            return Ok(());
        }
        if n == 1 {
            return Ok(());
        }
        let n = n as usize;
        let orig: Vec<ProteinSeq> = slf.inner.as_slice().to_vec();
        slf.inner.clear();
        slf.inner.reserve(orig.len() * n);
        for _ in 0..n {
            slf.inner.extend(orig.iter().cloned());
        }
        Ok(())
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DNABatch>()?;
    m.add_class::<RNABatch>()?;
    m.add_class::<ProteinBatch>()?;
    Ok(())
}
