#![allow(clippy::useless_conversion)]

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

use biorust_core::phylo;

use crate::msa::{AlignmentDNA, AlignmentProtein};

#[pyclass(frozen, name = "DistanceMatrix")]
pub struct PyDistanceMatrix {
    inner: phylo::DistanceMatrix,
}

#[pymethods]
impl PyDistanceMatrix {
    #[getter]
    fn n(&self) -> usize {
        self.inner.n()
    }

    fn labels(&self) -> Vec<String> {
        self.inner.labels().iter().map(|s| s.to_string()).collect()
    }

    fn get(&self, i: usize, j: usize) -> PyResult<f64> {
        let n = self.inner.n();
        if i >= n || j >= n {
            return Err(PyValueError::new_err(format!(
                "index ({}, {}) out of range for {}x{} matrix",
                i, j, n, n
            )));
        }
        Ok(self.inner.get(i, j))
    }

    fn to_list(&self) -> Vec<f64> {
        self.inner.data().to_vec()
    }

    fn to_list_of_lists(&self) -> Vec<Vec<f64>> {
        let n = self.inner.n();
        (0..n)
            .map(|i| (0..n).map(|j| self.inner.get(i, j)).collect())
            .collect()
    }

    fn __len__(&self) -> usize {
        self.inner.n()
    }

    fn __repr__(&self) -> String {
        format!("DistanceMatrix(n={})", self.inner.n())
    }

    fn __str__(&self) -> String {
        let n = self.inner.n();
        let labels = self.inner.labels();
        let pad = labels.iter().map(|l| l.len()).max().unwrap_or(0);

        let mut lines = Vec::with_capacity(n + 1);

        // Header line
        let mut header = format!("{:>pad$}", "", pad = pad + 2);
        for label in labels {
            header.push_str(&format!("{:>10}", &**label));
        }
        lines.push(header);

        for (i, label) in labels.iter().enumerate() {
            let mut row = format!("{:>pad$}  ", &**label, pad = pad);
            for j in 0..n {
                row.push_str(&format!("{:>10.4}", self.inner.get(i, j)));
            }
            lines.push(row);
        }

        lines.join("\n")
    }
}

#[pyclass(frozen, name = "PhyloTree")]
pub struct PyPhyloTree {
    inner: phylo::PhyloTree,
}

#[pymethods]
impl PyPhyloTree {
    fn to_newick(&self) -> String {
        phylo::to_newick(&self.inner)
    }

    fn ascii_diagram(&self) -> String {
        let start = self
            .inner
            .root()
            .unwrap_or_else(|| self.inner.num_nodes().saturating_sub(1));

        let mut out = String::new();
        out.push_str(&format_node_label(&self.inner, start));
        out.push('\n');

        let children = self.inner.node(start).children.clone();
        for (i, child) in children.iter().enumerate() {
            let last = i + 1 == children.len();
            write_ascii_subtree(&self.inner, *child, "", last, &mut out);
        }

        if out.ends_with('\n') {
            out.pop();
        }

        out
    }

    fn leaf_labels(&self) -> Vec<String> {
        self.inner.leaf_labels()
    }

    fn num_leaves(&self) -> usize {
        self.inner.num_leaves()
    }

    fn num_nodes(&self) -> usize {
        self.inner.num_nodes()
    }

    fn is_rooted(&self) -> bool {
        self.inner.root().is_some()
    }

    fn __repr__(&self) -> String {
        let kind = if self.inner.root().is_some() {
            "rooted"
        } else {
            "unrooted"
        };
        format!(
            "PhyloTree(leaves={}, nodes={}, {})",
            self.inner.num_leaves(),
            self.inner.num_nodes(),
            kind
        )
    }

    fn __str__(&self) -> String {
        self.to_newick()
    }
}

fn format_node_label(tree: &phylo::PhyloTree, idx: usize) -> String {
    let node = tree.node(idx);
    let mut label = if let Some(ref l) = node.label {
        l.to_string()
    } else if node.children.is_empty() {
        format!("leaf{}", idx)
    } else {
        format!("node{}", idx)
    };

    if let Some(bl) = node.branch_length {
        label.push_str(&format!(":{:.6}", bl));
    }

    label
}

fn write_ascii_subtree(
    tree: &phylo::PhyloTree,
    idx: usize,
    prefix: &str,
    is_last: bool,
    out: &mut String,
) {
    out.push_str(prefix);
    out.push_str(if is_last { "`-- " } else { "|-- " });
    out.push_str(&format_node_label(tree, idx));
    out.push('\n');

    let child_prefix = format!("{}{}", prefix, if is_last { "    " } else { "|   " });
    let children = tree.node(idx).children.clone();
    for (i, child) in children.iter().enumerate() {
        let last = i + 1 == children.len();
        write_ascii_subtree(tree, *child, &child_prefix, last, out);
    }
}

#[pyfunction]
#[pyo3(signature = (alignment, model = "p-distance"))]
fn distance_matrix(
    py: Python<'_>,
    alignment: &Bound<'_, PyAny>,
    model: &str,
) -> PyResult<PyDistanceMatrix> {
    // Try DNA first, then Protein
    if let Ok(dna) = alignment.extract::<PyRef<'_, AlignmentDNA>>() {
        let dna_model = match model {
            "p-distance" => phylo::DnaDistanceModel::PDistance,
            "jc69" => phylo::DnaDistanceModel::JukesCantor,
            "k2p" => phylo::DnaDistanceModel::Kimura2P,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "unknown DNA distance model '{}' (valid: 'p-distance', 'jc69', 'k2p')",
                    model
                )));
            }
        };

        let seqs_ref = dna.seqs_ref();
        let labels = dna.labels_cloned();
        let seq_bytes: Vec<&[u8]> = seqs_ref.iter().map(|s| s.as_bytes()).collect();

        let dm = py.allow_threads(|| phylo::dna_distance_matrix(&seq_bytes, labels, dna_model));

        return dm
            .map(|d| PyDistanceMatrix { inner: d })
            .map_err(|e| PyValueError::new_err(e.to_string()));
    }

    if let Ok(prot) = alignment.extract::<PyRef<'_, AlignmentProtein>>() {
        let prot_model = match model {
            "p-distance" => phylo::ProteinDistanceModel::PDistance,
            "poisson" => phylo::ProteinDistanceModel::Poisson,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "unknown protein distance model '{}' (valid: 'p-distance', 'poisson')",
                    model
                )));
            }
        };

        let seqs_ref = prot.seqs_ref();
        let labels = prot.labels_cloned();
        let seq_bytes: Vec<&[u8]> = seqs_ref.iter().map(|s| s.as_bytes()).collect();

        let dm =
            py.allow_threads(|| phylo::protein_distance_matrix(&seq_bytes, labels, prot_model));

        return dm
            .map(|d| PyDistanceMatrix { inner: d })
            .map_err(|e| PyValueError::new_err(e.to_string()));
    }

    Err(PyValueError::new_err(
        "alignment must be AlignmentDNA or AlignmentProtein",
    ))
}

#[pyfunction]
#[pyo3(signature = (dist_matrix, method = "nj"))]
fn build_tree(
    py: Python<'_>,
    dist_matrix: &PyDistanceMatrix,
    method: &str,
) -> PyResult<PyPhyloTree> {
    let dm = &dist_matrix.inner;

    let tree = match method {
        "nj" => py.allow_threads(|| phylo::neighbor_joining(dm)),
        "upgma" => py.allow_threads(|| phylo::upgma(dm)),
        _ => {
            return Err(PyValueError::new_err(format!(
                "unknown tree method '{}' (valid: 'nj', 'upgma')",
                method
            )));
        }
    };

    tree.map(|t| PyPhyloTree { inner: t })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDistanceMatrix>()?;
    m.add_class::<PyPhyloTree>()?;
    m.add_function(wrap_pyfunction!(distance_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(build_tree, m)?)?;
    Ok(())
}
