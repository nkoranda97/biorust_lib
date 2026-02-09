pub mod encode;
pub mod global_simd;
pub mod local_simd;
pub mod matrices;
pub mod scalar_ref;
mod simd_utils;
pub mod types;

pub use encode::{encode_dna, encode_protein, EncodedSeq};
pub use types::{AlignmentMode, AlignmentResult, Cigar, CigarOp, Scoring};

#[cfg(test)]
mod tests;

/// Conservative limit to prevent i16 overflow in SIMD kernels.
/// SIMD implementations use i16 for DP values, so we need to ensure
/// max_score * sequence_length stays well below i16::MAX (32767).
const SIMD_MAX_SAFE_SCORE: i32 = 30000;

fn simd_safe_len(len: usize, scoring: &Scoring) -> bool {
    let max_abs = scoring.max_abs_score();
    let bound = max_abs.saturating_mul(len as i32);
    bound <= SIMD_MAX_SAFE_SCORE
}

#[allow(dead_code)]
pub(crate) fn score_alignment_from_cigar(
    query: &[u8],
    target: &[u8],
    cigar: &Cigar,
    scoring: &Scoring,
) -> f32 {
    let mut qi = 0usize;
    let mut ti = 0usize;
    let mut score = 0.0f32;
    for (op, len) in &cigar.ops {
        match op {
            CigarOp::Match => {
                for _ in 0..*len {
                    score += scoring.score(query[qi], target[ti]) as f32;
                    qi += 1;
                    ti += 1;
                }
            }
            CigarOp::Ins => {
                let is_terminal = scoring.end_gap
                    && ((qi == 0 && ti == 0) || (ti == target.len() && qi + len == query.len()));
                let (gap_open, gap_extend) = if is_terminal {
                    (scoring.end_gap_open, scoring.end_gap_extend)
                } else {
                    (scoring.gap_open, scoring.gap_extend)
                };
                if *len > 0 {
                    score += gap_open + gap_extend * (*len as f32 - 1.0);
                }
                qi += *len;
            }
            CigarOp::Del => {
                let is_terminal = scoring.end_gap
                    && ((qi == 0 && ti == 0) || (qi == query.len() && ti + len == target.len()));
                let (gap_open, gap_extend) = if is_terminal {
                    (scoring.end_gap_open, scoring.end_gap_extend)
                } else {
                    (scoring.gap_open, scoring.gap_extend)
                };
                if *len > 0 {
                    score += gap_open + gap_extend * (*len as f32 - 1.0);
                }
                ti += *len;
            }
        }
    }
    score
}

pub fn align_local(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
    traceback: bool,
) -> AlignmentResult {
    if scoring.matrix.is_some() {
        assert_eq!(
            scoring
                .alphabet_size
                .expect("alphabet_size must be set when matrix is present"),
            query.alphabet_size,
            "scoring matrix alphabet size mismatch"
        );
    }

    if !traceback
        && scoring.simd_compatible()
        && simd_safe_len(query.len().max(target.len()), scoring)
    {
        #[cfg(feature = "simd")]
        {
            let (score, end_q, end_t) = local_simd::align_local_score(query, target, scoring);
            return AlignmentResult {
                score,
                query_end: end_q,
                target_end: end_t,
                query_start: None,
                target_start: None,
                cigar: None,
            };
        }
    }

    scalar_ref::align_local_scalar(query, target, scoring, traceback)
}

pub fn align_global(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
    traceback: bool,
) -> AlignmentResult {
    if scoring.matrix.is_some() {
        assert_eq!(
            scoring
                .alphabet_size
                .expect("alphabet_size must be set when matrix is present"),
            query.alphabet_size,
            "scoring matrix alphabet size mismatch"
        );
    }

    if !traceback
        && scoring.simd_compatible()
        && simd_safe_len(query.len().max(target.len()), scoring)
    {
        #[cfg(feature = "simd")]
        {
            let (score, end_q, end_t) = global_simd::align_global_score(query, target, scoring);
            return AlignmentResult {
                score,
                query_end: end_q,
                target_end: end_t,
                query_start: Some(0),
                target_start: Some(0),
                cigar: None,
            };
        }
    }

    scalar_ref::align_global_scalar(query, target, scoring, traceback)
}
