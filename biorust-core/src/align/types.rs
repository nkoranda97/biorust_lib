#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignmentMode {
    Local,
    Global,
}

/// CIGAR operations consume sequence coordinates.
/// Ins consumes query (gap in target), Del consumes target (gap in query).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    /// Consumes query and target.
    Match,
    /// Consumes query, gap in target.
    Ins,
    /// Consumes target, gap in query.
    Del,
}

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Cigar {
    pub(crate) ops: Vec<(CigarOp, usize)>,
}

impl Cigar {
    pub fn ops(&self) -> &[(CigarOp, usize)] {
        &self.ops
    }

    pub fn into_ops(self) -> Vec<(CigarOp, usize)> {
        self.ops
    }

    pub fn push(&mut self, op: CigarOp, len: usize) {
        if len == 0 {
            return;
        }
        if let Some((last_op, last_len)) = self.ops.last_mut() {
            if *last_op == op {
                *last_len += len;
                return;
            }
        }
        self.ops.push((op, len));
    }

    pub fn len(&self) -> usize {
        self.ops.iter().map(|(_, n)| *n).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }
}

impl std::fmt::Display for Cigar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for &(op, len) in &self.ops {
            let ch = match op {
                CigarOp::Match => 'M',
                CigarOp::Ins => 'I',
                CigarOp::Del => 'D',
            };
            write!(f, "{len}{ch}")?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentResult {
    pub score: f32,
    pub query_end: usize,
    pub target_end: usize,
    pub query_start: Option<usize>,
    pub target_start: Option<usize>,
    pub cigar: Option<Cigar>,
}

#[derive(Clone, Debug)]
pub struct Scoring {
    pub(crate) match_score: i16,
    pub(crate) mismatch_score: i16,
    pub(crate) gap_open: f32,
    pub(crate) gap_extend: f32,
    pub(crate) matrix: Option<Vec<i16>>,
    pub(crate) alphabet_size: Option<usize>,
    pub(crate) end_gap: bool,
    pub(crate) end_gap_open: f32,
    pub(crate) end_gap_extend: f32,
}

impl Scoring {
    pub fn gap_open(&self) -> f32 {
        self.gap_open
    }

    pub fn gap_extend(&self) -> f32 {
        self.gap_extend
    }

    pub fn end_gap(&self) -> bool {
        self.end_gap
    }

    pub fn end_gap_open(&self) -> f32 {
        self.end_gap_open
    }

    pub fn end_gap_extend(&self) -> f32 {
        self.end_gap_extend
    }

    pub fn match_score(&self) -> i16 {
        self.match_score
    }

    pub fn mismatch_score(&self) -> i16 {
        self.mismatch_score
    }

    pub fn matrix(&self) -> Option<&[i16]> {
        self.matrix.as_deref()
    }

    pub fn alphabet_size_opt(&self) -> Option<usize> {
        self.alphabet_size
    }

    pub fn simple(
        match_score: i16,
        mismatch_score: i16,
        gap_open: f32,
        gap_extend: f32,
    ) -> crate::error::BioResult<Self> {
        if gap_open > 0.0 {
            return Err(crate::error::BioError::InvalidScoring {
                msg: format!("gap_open must be <= 0, got {gap_open}"),
            });
        }
        if gap_extend > 0.0 {
            return Err(crate::error::BioError::InvalidScoring {
                msg: format!("gap_extend must be <= 0, got {gap_extend}"),
            });
        }
        Ok(Self {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            matrix: None,
            alphabet_size: None,
            end_gap: false,
            end_gap_open: gap_open,
            end_gap_extend: gap_extend,
        })
    }

    pub fn with_matrix(
        matrix: Vec<i16>,
        alphabet_size: usize,
        gap_open: f32,
        gap_extend: f32,
    ) -> crate::error::BioResult<Self> {
        if gap_open > 0.0 {
            return Err(crate::error::BioError::InvalidScoring {
                msg: format!("gap_open must be <= 0, got {gap_open}"),
            });
        }
        if gap_extend > 0.0 {
            return Err(crate::error::BioError::InvalidScoring {
                msg: format!("gap_extend must be <= 0, got {gap_extend}"),
            });
        }
        if alphabet_size == 0 {
            return Err(crate::error::BioError::InvalidScoring {
                msg: "alphabet_size must be > 0".into(),
            });
        }
        if matrix.len() != alphabet_size * alphabet_size {
            return Err(crate::error::BioError::InvalidScoring {
                msg: format!(
                    "matrix length {} doesn't match alphabet_size² {}",
                    matrix.len(),
                    alphabet_size * alphabet_size
                ),
            });
        }
        Ok(Self {
            match_score: 0,
            mismatch_score: 0,
            gap_open,
            gap_extend,
            matrix: Some(matrix),
            alphabet_size: Some(alphabet_size),
            end_gap: false,
            end_gap_open: gap_open,
            end_gap_extend: gap_extend,
        })
    }

    pub fn with_end_gaps(mut self, end_gap_open: f32, end_gap_extend: f32) -> Self {
        self.end_gap = true;
        self.end_gap_open = end_gap_open;
        self.end_gap_extend = end_gap_extend;
        self
    }

    pub fn simd_compatible(&self) -> bool {
        if self.end_gap {
            return false;
        }
        let gap_open = self.gap_open;
        let gap_extend = self.gap_extend;
        if !gap_open.is_finite() || !gap_extend.is_finite() {
            return false;
        }
        let gap_open_int = gap_open.fract().abs() < f32::EPSILON;
        let gap_extend_int = gap_extend.fract().abs() < f32::EPSILON;
        if !gap_open_int || !gap_extend_int {
            return false;
        }
        if gap_open > 0.0 || gap_extend > 0.0 {
            return false;
        }
        let max_i16 = i16::MAX as f32;
        gap_open.abs() <= max_i16 && gap_extend.abs() <= max_i16
    }

    pub fn gap_open_i16(&self) -> i16 {
        (-self.gap_open).round() as i16
    }

    pub fn gap_extend_i16(&self) -> i16 {
        (-self.gap_extend).round() as i16
    }

    #[inline]
    pub fn score(&self, a: u8, b: u8) -> i16 {
        if let Some(matrix) = &self.matrix {
            let alpha = self
                .alphabet_size
                .expect("alphabet_size must be set when matrix is present");
            let idx = (a as usize) * alpha + (b as usize);
            matrix[idx]
        } else if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }

    pub fn max_abs_score(&self) -> i32 {
        let mut max_abs = self.match_score.abs().max(self.mismatch_score.abs()) as i32;
        if let Some(matrix) = &self.matrix {
            for &v in matrix {
                max_abs = max_abs.max(v.abs() as i32);
            }
        }
        max_abs = max_abs.max(self.gap_open.abs().ceil() as i32);
        max_abs = max_abs.max(self.gap_extend.abs().ceil() as i32);
        if self.end_gap {
            max_abs = max_abs.max(self.end_gap_open.abs().ceil() as i32);
            max_abs = max_abs.max(self.end_gap_extend.abs().ceil() as i32);
        }
        max_abs
    }
}
