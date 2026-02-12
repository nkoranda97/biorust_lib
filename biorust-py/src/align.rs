#![allow(clippy::useless_conversion)]

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use biorust_core::align as core_align;

use crate::dna::DNA;
use crate::protein::Protein;

#[pyclass(frozen)]
pub struct Scoring {
    pub(crate) inner: core_align::Scoring,
    use_matrix: bool,
}

fn to_i16(name: &str, value: i64) -> PyResult<i16> {
    i16::try_from(value).map_err(|_| {
        PyValueError::new_err(format!("{name} must fit in int16 range (-32768..=32767)"))
    })
}

fn to_f32(name: &str, value: f64) -> PyResult<f32> {
    let v = value as f32;
    if !v.is_finite() {
        return Err(PyValueError::new_err(format!(
            "{name} must be a finite float"
        )));
    }
    Ok(v)
}

#[pymethods]
impl Scoring {
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (match_score=None, mismatch_score=None, gap_open=-2.0, gap_extend=-1.0, matrix=None, alphabet_size=None, end_gap=true, end_gap_open=None, end_gap_extend=None, use_matrix=None))]
    fn new(
        match_score: Option<i64>,
        mismatch_score: Option<i64>,
        gap_open: f64,
        gap_extend: f64,
        matrix: Option<&Bound<'_, PyAny>>,
        alphabet_size: Option<usize>,
        end_gap: bool,
        end_gap_open: Option<f64>,
        end_gap_extend: Option<f64>,
        use_matrix: Option<bool>,
    ) -> PyResult<Self> {
        let gap_open = to_f32("gap_open", gap_open)?;
        let gap_extend = to_f32("gap_extend", gap_extend)?;
        let end_gap_open = match end_gap_open {
            Some(v) => to_f32("end_gap_open", v)?,
            None => 0.0,
        };
        let end_gap_extend = match end_gap_extend {
            Some(v) => to_f32("end_gap_extend", v)?,
            None => 0.0,
        };
        // If the user explicitly provided match/mismatch scores, use simple scoring.
        // If neither was provided and no explicit matrix, auto-select EDNAFULL/BLOSUM62.
        let scores_provided = match_score.is_some() || mismatch_score.is_some();
        let match_score = match_score.unwrap_or(2);
        let mismatch_score = mismatch_score.unwrap_or(-1);

        let mut scoring = if let Some(matrix) = matrix {
            if let Ok(name) = matrix.extract::<String>() {
                let def = core_align::matrices::matrix_by_name(&name)
                    .ok_or_else(|| PyValueError::new_err(format!("unknown matrix name: {name}")))?;
                core_align::Scoring::with_matrix(
                    def.scores.to_vec(),
                    def.alphabet.len(),
                    gap_open,
                    gap_extend,
                )
                .map_err(|e| PyValueError::new_err(e.to_string()))?
            } else if let Ok(matrix_vals) = matrix.extract::<Vec<i64>>() {
                let alpha = alphabet_size.ok_or_else(|| {
                    PyValueError::new_err("alphabet_size is required when matrix is provided")
                })?;
                if alpha == 0 {
                    return Err(PyValueError::new_err("alphabet_size must be > 0"));
                }
                if matrix_vals.len() != alpha * alpha {
                    return Err(PyValueError::new_err(format!(
                        "matrix length must be alphabet_size^2 (got {}, expected {})",
                        matrix_vals.len(),
                        alpha * alpha
                    )));
                }
                let mut mtx = Vec::with_capacity(matrix_vals.len());
                for (i, v) in matrix_vals.into_iter().enumerate() {
                    let val = to_i16("matrix value", v).map_err(|_| {
                        PyValueError::new_err(format!(
                            "matrix value at index {i} must fit in int16 range"
                        ))
                    })?;
                    mtx.push(val);
                }
                core_align::Scoring::with_matrix(mtx, alpha, gap_open, gap_extend)
                    .map_err(|e| PyValueError::new_err(e.to_string()))?
            } else {
                return Err(PyTypeError::new_err(
                    "matrix must be a list of ints or a known matrix name",
                ));
            }
        } else {
            core_align::Scoring::simple(
                to_i16("match_score", match_score)?,
                to_i16("mismatch_score", mismatch_score)?,
                gap_open,
                gap_extend,
            )
            .map_err(|e| PyValueError::new_err(e.to_string()))?
        };

        if !end_gap {
            scoring = scoring
                .with_end_gaps(end_gap_open, end_gap_extend)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        }

        // use_matrix resolution:
        // - Explicit use_matrix=True/False from user: honor it
        // - matrix= provided: always use it
        // - match_score/mismatch_score provided: use simple scoring
        // - Nothing provided: auto-select EDNAFULL/BLOSUM62
        let use_matrix = match use_matrix {
            Some(v) => v || matrix.is_some(),
            None => matrix.is_some() || !scores_provided,
        };

        Ok(Self {
            inner: scoring,
            use_matrix,
        })
    }

    #[staticmethod]
    #[allow(clippy::useless_conversion)]
    #[pyo3(signature = (matrix, alphabet_size, gap_open=-2.0, gap_extend=-1.0, end_gap=true, end_gap_open=None, end_gap_extend=None))]
    fn with_matrix(
        matrix: Vec<i64>,
        alphabet_size: usize,
        gap_open: f64,
        gap_extend: f64,
        end_gap: bool,
        end_gap_open: Option<f64>,
        end_gap_extend: Option<f64>,
    ) -> PyResult<Self> {
        if alphabet_size == 0 {
            return Err(PyValueError::new_err("alphabet_size must be > 0"));
        }
        if matrix.len() != alphabet_size * alphabet_size {
            return Err(PyValueError::new_err(format!(
                "matrix length must be alphabet_size^2 (got {}, expected {})",
                matrix.len(),
                alphabet_size * alphabet_size
            )));
        }
        let mut mtx = Vec::with_capacity(matrix.len());
        for (i, v) in matrix.into_iter().enumerate() {
            let val = to_i16("matrix value", v).map_err(|_| {
                PyValueError::new_err(format!("matrix value at index {i} must fit in int16 range"))
            })?;
            mtx.push(val);
        }
        let gap_open = to_f32("gap_open", gap_open)?;
        let gap_extend = to_f32("gap_extend", gap_extend)?;
        let end_gap_open = match end_gap_open {
            Some(v) => to_f32("end_gap_open", v)?,
            None => 0.0,
        };
        let end_gap_extend = match end_gap_extend {
            Some(v) => to_f32("end_gap_extend", v)?,
            None => 0.0,
        };
        let mut scoring =
            core_align::Scoring::with_matrix(mtx, alphabet_size, gap_open, gap_extend)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        if !end_gap {
            scoring = scoring
                .with_end_gaps(end_gap_open, end_gap_extend)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        }
        Ok(Self {
            inner: scoring,
            use_matrix: true,
        })
    }

    #[staticmethod]
    fn matrix_names() -> Vec<&'static str> {
        core_align::matrices::matrix_names().to_vec()
    }

    fn __repr__(&self) -> String {
        if self.inner.matrix().is_some() {
            if self.inner.end_gap() {
                format!(
                    "Scoring(matrix=..., alphabet_size={}, gap_open={}, gap_extend={}, end_gap=false, end_gap_open={}, end_gap_extend={})",
                    self.inner
                        .alphabet_size_opt()
                        .expect("alphabet_size must be set when matrix is present"),
                    self.inner.gap_open(),
                    self.inner.gap_extend(),
                    self.inner.end_gap_open(),
                    self.inner.end_gap_extend()
                )
            } else {
                format!(
                    "Scoring(matrix=..., alphabet_size={}, gap_open={}, gap_extend={})",
                    self.inner
                        .alphabet_size_opt()
                        .expect("alphabet_size must be set when matrix is present"),
                    self.inner.gap_open(),
                    self.inner.gap_extend()
                )
            }
        } else if self.inner.end_gap() {
            format!(
                "Scoring(match_score={}, mismatch_score={}, gap_open={}, gap_extend={}, end_gap=false, end_gap_open={}, end_gap_extend={}, use_matrix={})",
                self.inner.match_score(),
                self.inner.mismatch_score(),
                self.inner.gap_open(),
                self.inner.gap_extend(),
                self.inner.end_gap_open(),
                self.inner.end_gap_extend(),
                self.use_matrix
            )
        } else {
            format!(
                "Scoring(match_score={}, mismatch_score={}, gap_open={}, gap_extend={}, use_matrix={})",
                self.inner.match_score(),
                self.inner.mismatch_score(),
                self.inner.gap_open(),
                self.inner.gap_extend(),
                self.use_matrix
            )
        }
    }
}

#[pyclass(frozen)]
pub struct AlignmentResult {
    inner: core_align::AlignmentResult,
    query: Vec<u8>,
    target: Vec<u8>,
}

fn cigar_to_py(cigar: &core_align::Cigar) -> Vec<(String, usize)> {
    cigar
        .ops()
        .iter()
        .map(|(op, len)| {
            let code = match op {
                core_align::CigarOp::Match => "M",
                core_align::CigarOp::Ins => "I",
                core_align::CigarOp::Del => "D",
            };
            (code.to_string(), *len)
        })
        .collect()
}

#[pymethods]
impl AlignmentResult {
    #[getter]
    fn score(&self) -> f32 {
        self.inner.score
    }

    #[getter]
    fn query_end(&self) -> usize {
        self.inner.query_end
    }

    #[getter]
    fn target_end(&self) -> usize {
        self.inner.target_end
    }

    #[getter]
    fn query_start(&self) -> Option<usize> {
        self.inner.query_start
    }

    #[getter]
    fn target_start(&self) -> Option<usize> {
        self.inner.target_start
    }

    #[getter]
    fn cigar(&self) -> Option<Vec<(String, usize)>> {
        self.inner.cigar.as_ref().map(cigar_to_py)
    }

    fn aligned_strings(&self) -> PyResult<(String, String)> {
        let (q_out, _mid_out, t_out) = self.alignment_parts()?;
        Ok((q_out, t_out))
    }

    fn alignment_diagram(&self) -> PyResult<String> {
        let (q_out, mid_out, t_out) = self.alignment_parts()?;
        Ok(format!("{q_out}\n{mid_out}\n{t_out}"))
    }

    fn __repr__(&self) -> String {
        let cigar = self
            .inner
            .cigar
            .as_ref()
            .map(|c| c.ops().len())
            .unwrap_or(0);
        format!(
            "AlignmentResult(score={}, query_end={}, target_end={}, cigar_ops={})",
            self.inner.score, self.inner.query_end, self.inner.target_end, cigar
        )
    }
}

impl AlignmentResult {
    fn alignment_parts(&self) -> PyResult<(String, String, String)> {
        let cigar = self
            .inner
            .cigar
            .as_ref()
            .ok_or_else(|| PyValueError::new_err("traceback was not requested"))?;
        let q_start = self
            .inner
            .query_start
            .ok_or_else(|| PyValueError::new_err("query_start is missing"))?;
        let t_start = self
            .inner
            .target_start
            .ok_or_else(|| PyValueError::new_err("target_start is missing"))?;

        let q_bytes = &self.query;
        let t_bytes = &self.target;

        let mut q_idx = q_start;
        let mut t_idx = t_start;
        let mut q_out = String::new();
        let mut mid_out = String::new();
        let mut t_out = String::new();

        for (op, len) in cigar.ops() {
            match op {
                core_align::CigarOp::Match => {
                    // Use checked_add to prevent overflow before bounds check
                    let q_end = q_idx
                        .checked_add(*len)
                        .ok_or_else(|| PyValueError::new_err("cigar length overflow"))?;
                    let t_end = t_idx
                        .checked_add(*len)
                        .ok_or_else(|| PyValueError::new_err("cigar length overflow"))?;
                    if q_end > q_bytes.len() || t_end > t_bytes.len() {
                        return Err(PyValueError::new_err("cigar exceeds sequence bounds"));
                    }
                    for k in 0..*len {
                        let qb = q_bytes[q_idx + k];
                        let tb = t_bytes[t_idx + k];
                        q_out.push(qb as char);
                        t_out.push(tb as char);
                        mid_out.push(if qb == tb { '|' } else { '*' });
                    }
                    q_idx += len;
                    t_idx += len;
                }
                core_align::CigarOp::Ins => {
                    let q_end = q_idx
                        .checked_add(*len)
                        .ok_or_else(|| PyValueError::new_err("cigar length overflow"))?;
                    if q_end > q_bytes.len() {
                        return Err(PyValueError::new_err("cigar exceeds sequence bounds"));
                    }
                    for k in 0..*len {
                        let qb = q_bytes[q_idx + k];
                        q_out.push(qb as char);
                        t_out.push('-');
                        mid_out.push(' ');
                    }
                    q_idx += len;
                }
                core_align::CigarOp::Del => {
                    let t_end = t_idx
                        .checked_add(*len)
                        .ok_or_else(|| PyValueError::new_err("cigar length overflow"))?;
                    if t_end > t_bytes.len() {
                        return Err(PyValueError::new_err("cigar exceeds sequence bounds"));
                    }
                    for k in 0..*len {
                        let tb = t_bytes[t_idx + k];
                        q_out.push('-');
                        t_out.push(tb as char);
                        mid_out.push(' ');
                    }
                    t_idx += len;
                }
            }
        }

        Ok((q_out, mid_out, t_out))
    }
}

enum SeqKind {
    Dna(Vec<u8>),
    Protein(Vec<u8>),
}

fn extract_seq(obj: &Bound<'_, PyAny>) -> PyResult<SeqKind> {
    if let Ok(dna) = obj.extract::<PyRef<'_, DNA>>() {
        return Ok(SeqKind::Dna(dna.as_bytes().to_vec()));
    }
    if let Ok(protein) = obj.extract::<PyRef<'_, Protein>>() {
        return Ok(SeqKind::Protein(protein.as_bytes().to_vec()));
    }
    Err(PyTypeError::new_err(
        "query/target must be DNA or Protein objects",
    ))
}

fn align_internal(
    query: &Bound<'_, PyAny>,
    target: &Bound<'_, PyAny>,
    scoring: &Scoring,
    traceback: bool,
    local: bool,
) -> PyResult<AlignmentResult> {
    let q = extract_seq(query)?;
    let t = extract_seq(target)?;

    let (q_enc, t_enc, is_dna, q_bytes, t_bytes) = match (q, t) {
        (SeqKind::Dna(q), SeqKind::Dna(t)) => (
            core_align::encode_dna(&q).map_err(|e| PyValueError::new_err(e.to_string()))?,
            core_align::encode_dna(&t).map_err(|e| PyValueError::new_err(e.to_string()))?,
            true,
            q,
            t,
        ),
        (SeqKind::Protein(q), SeqKind::Protein(t)) => (
            core_align::encode_protein(&q).map_err(|e| PyValueError::new_err(e.to_string()))?,
            core_align::encode_protein(&t).map_err(|e| PyValueError::new_err(e.to_string()))?,
            false,
            q,
            t,
        ),
        _ => {
            return Err(PyValueError::new_err(
                "query and target must be the same sequence type",
            ))
        }
    };

    let auto_scoring = if scoring.inner.matrix().is_none() && scoring.use_matrix {
        let def = if is_dna {
            core_align::matrices::matrix_by_name("EDNAFULL").expect("EDNAFULL matrix is available")
        } else {
            core_align::matrices::matrix_by_name("BLOSUM62").expect("BLOSUM62 matrix is available")
        };
        if def.alphabet.len() != q_enc.alphabet_size() {
            return Err(PyValueError::new_err(
                "scoring matrix alphabet size does not match sequence alphabet",
            ));
        }
        let mut sc = core_align::Scoring::with_matrix(
            def.scores.to_vec(),
            def.alphabet.len(),
            scoring.inner.gap_open(),
            scoring.inner.gap_extend(),
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        if scoring.inner.end_gap() {
            sc = sc
                .with_end_gaps(scoring.inner.end_gap_open(), scoring.inner.end_gap_extend())
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        }
        Some(sc)
    } else {
        None
    };
    let scoring_ref = auto_scoring.as_ref().unwrap_or(&scoring.inner);
    if scoring_ref.matrix().is_some()
        && scoring_ref
            .alphabet_size_opt()
            .expect("alphabet_size must be set when matrix is present")
            != q_enc.alphabet_size()
    {
        return Err(PyValueError::new_err(
            "scoring matrix alphabet size does not match sequence alphabet",
        ));
    }

    // Release GIL during alignment computation to allow other Python threads to run
    let py = query.py();
    let inner = py.allow_threads(|| {
        if local {
            core_align::align_local(&q_enc, &t_enc, scoring_ref, traceback)
        } else {
            core_align::align_global(&q_enc, &t_enc, scoring_ref, traceback)
        }
    });

    Ok(AlignmentResult {
        inner,
        query: q_bytes,
        target: t_bytes,
    })
}

#[pyfunction]
#[pyo3(signature = (query, target, scoring, traceback=false))]
#[allow(clippy::useless_conversion)]
fn align_local(
    query: &Bound<'_, PyAny>,
    target: &Bound<'_, PyAny>,
    scoring: &Scoring,
    traceback: bool,
) -> PyResult<AlignmentResult> {
    align_internal(query, target, scoring, traceback, true)
}

#[pyfunction]
#[pyo3(signature = (query, target, scoring, traceback=false))]
#[allow(clippy::useless_conversion)]
fn align_global(
    query: &Bound<'_, PyAny>,
    target: &Bound<'_, PyAny>,
    scoring: &Scoring,
    traceback: bool,
) -> PyResult<AlignmentResult> {
    align_internal(query, target, scoring, traceback, false)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Scoring>()?;
    m.add_class::<AlignmentResult>()?;
    m.add_function(wrap_pyfunction!(align_local, m)?)?;
    m.add_function(wrap_pyfunction!(align_global, m)?)?;
    Ok(())
}
