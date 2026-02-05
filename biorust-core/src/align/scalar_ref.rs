//! Alignment DP uses i for target (rows) and j for query (columns).

use super::encode::EncodedSeq;
use super::types::{AlignmentResult, Cigar, CigarOp, Scoring};

// Bits 0-1: H direction
const DIR_DIAG: u8 = 0;
const DIR_INS: u8 = 1; // F (horizontal gap, query insertion)
const DIR_DEL: u8 = 2; // E (vertical gap, target deletion)
const DIR_ZERO: u8 = 3;

// Traceback for E/F: 0 = opened from H, 1 = extended from E/F.
const TRACE_E_FROM_H: u8 = 0;
const TRACE_E_FROM_E: u8 = 1;
const TRACE_F_FROM_H: u8 = 0;
const TRACE_F_FROM_F: u8 = 1;

fn push_rev(ops: &mut Vec<(CigarOp, usize)>, op: CigarOp, len: usize) {
    if len == 0 {
        return;
    }
    if let Some((last_op, last_len)) = ops.last_mut() {
        if *last_op == op {
            *last_len += len;
            return;
        }
    }
    ops.push((op, len));
}

fn finalize_cigar(rev_ops: Vec<(CigarOp, usize)>) -> Cigar {
    // Reverse without merging — the traceback already encoded gap block boundaries
    Cigar {
        ops: rev_ops.into_iter().rev().collect(),
    }
}

pub fn align_local_scalar(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
    traceback: bool,
) -> AlignmentResult {
    let m = query.codes.len();
    let n = target.codes.len();
    if m == 0 || n == 0 {
        return AlignmentResult {
            score: 0.0,
            query_end: 0,
            target_end: 0,
            query_start: Some(0),
            target_start: Some(0),
            cigar: Some(Cigar::default()),
        };
    }

    let neg_inf: f32 = f32::NEG_INFINITY;
    let gap_open = scoring.gap_open;
    let gap_extend = scoring.gap_extend;

    let mut h_row = vec![0f32; m + 1];
    let mut e_row = vec![neg_inf; m + 1];
    let mut trace_h = if traceback {
        vec![DIR_ZERO; (n + 1) * (m + 1)]
    } else {
        Vec::new()
    };
    let mut trace_e = if traceback {
        vec![TRACE_E_FROM_H; (n + 1) * (m + 1)]
    } else {
        Vec::new()
    };
    let mut trace_f = if traceback {
        vec![TRACE_F_FROM_H; (n + 1) * (m + 1)]
    } else {
        Vec::new()
    };

    let mut max_score = 0f32;
    let mut end_i = 0usize;
    let mut end_j = 0usize;

    for i in 1..=n {
        let t = target.codes[i - 1];
        let mut h_diag = 0f32;
        let mut f = neg_inf;
        h_row[0] = 0.0;
        if traceback {
            trace_h[i * (m + 1)] = DIR_ZERO;
        }
        for j in 1..=m {
            let h_up = h_row[j];
            let e_open = h_up + gap_open;
            let e_ext = e_row[j] + gap_extend;
            let e_from_ext = e_ext > e_open;
            e_row[j] = if e_from_ext { e_ext } else { e_open };
            let f_open = h_row[j - 1] + gap_open;
            let f_ext = f + gap_extend;
            let f_from_ext = f_ext > f_open;
            f = if f_from_ext { f_ext } else { f_open };
            let score_diag = h_diag + scoring.score(query.codes[j - 1], t) as f32;
            let mut h = score_diag;
            let mut d = DIR_DIAG;
            // Tie-breaking policy:
            // DIAG > DEL > INS (because we use strict > comparisons)
            // Multiple optimal alignments may exist; this is intentional.
            if e_row[j] > h {
                h = e_row[j];
                d = DIR_DEL;
            }
            if f > h {
                h = f;
                d = DIR_INS;
            }
            if h < 0.0 {
                h = 0.0;
                d = DIR_ZERO;
            }
            if traceback {
                let idx = i * (m + 1) + j;
                // trace_h encodes the predecessor of H(i,j); trace_e/f encode the predecessor of E/F.
                trace_e[idx] = if e_from_ext {
                    TRACE_E_FROM_E
                } else {
                    TRACE_E_FROM_H
                };
                trace_f[idx] = if f_from_ext {
                    TRACE_F_FROM_F
                } else {
                    TRACE_F_FROM_H
                };
                trace_h[idx] = d;
            }
            h_row[j] = h;
            if h > max_score {
                max_score = h;
                end_i = i;
                end_j = j;
            }
            h_diag = h_up;
        }
    }

    if !traceback {
        return AlignmentResult {
            score: max_score,
            query_end: end_j.saturating_sub(1),
            target_end: end_i.saturating_sub(1),
            query_start: None,
            target_start: None,
            cigar: None,
        };
    }

    let mut i = end_i;
    let mut j = end_j;
    let mut rev_ops: Vec<(CigarOp, usize)> = Vec::new();
    // State: 0 = in H matrix, 1 = in E matrix (vertical/DEL), 2 = in F matrix (horizontal/INS)
    let mut state = 0u8;

    loop {
        if i == 0 && j == 0 {
            break;
        }
        match state {
            0 => {
                let d = trace_h[i * (m + 1) + j];
                match d {
                    DIR_ZERO => break,
                    DIR_DIAG => {
                        push_rev(&mut rev_ops, CigarOp::Match, 1);
                        i -= 1;
                        j -= 1;
                    }
                    DIR_DEL => {
                        state = 1;
                    }
                    DIR_INS => {
                        state = 2;
                    }
                    _ => break,
                }
            }
            1 => {
                if i == 0 {
                    break;
                }
                let d = trace_e[i * (m + 1) + j];
                let extending = d == TRACE_E_FROM_E;
                push_rev(&mut rev_ops, CigarOp::Del, 1);
                i -= 1;
                if !extending {
                    state = 0;
                }
            }
            2 => {
                if j == 0 {
                    break;
                }
                let d = trace_f[i * (m + 1) + j];
                let extending = d == TRACE_F_FROM_F;
                push_rev(&mut rev_ops, CigarOp::Ins, 1);
                j -= 1;
                if !extending {
                    state = 0;
                }
            }
            _ => break,
        }
    }

    AlignmentResult {
        score: max_score,
        query_end: end_j.saturating_sub(1),
        target_end: end_i.saturating_sub(1),
        query_start: Some(j),
        target_start: Some(i),
        cigar: Some(finalize_cigar(rev_ops)),
    }
}

pub fn align_global_scalar(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
    traceback: bool,
) -> AlignmentResult {
    let m = query.codes.len();
    let n = target.codes.len();
    let neg_inf: f32 = f32::NEG_INFINITY;
    let gap_open = scoring.gap_open;
    let gap_extend = scoring.gap_extend;
    let end_gap_open = if scoring.end_gap {
        scoring.end_gap_open
    } else {
        gap_open
    };
    let end_gap_extend = if scoring.end_gap {
        scoring.end_gap_extend
    } else {
        gap_extend
    };

    if m == 0 || n == 0 {
        let len = if m == 0 { n } else { m };
        let score = if len == 0 {
            0.0
        } else {
            end_gap_open + end_gap_extend * (len as f32 - 1.0)
        };
        let mut cigar = Cigar::default();
        if len > 0 {
            if m == 0 {
                cigar.push(CigarOp::Del, len);
            } else {
                cigar.push(CigarOp::Ins, len);
            }
        }
        return AlignmentResult {
            score,
            query_end: m.saturating_sub(1),
            target_end: n.saturating_sub(1),
            query_start: Some(0),
            target_start: Some(0),
            cigar: Some(cigar),
        };
    }

    let mut h_row = vec![0f32; m + 1];
    let mut e_row = vec![neg_inf; m + 1];
    let mut trace_h = if traceback {
        vec![DIR_DIAG; (n + 1) * (m + 1)]
    } else {
        Vec::new()
    };
    let mut trace_e = if traceback {
        vec![TRACE_E_FROM_H; (n + 1) * (m + 1)]
    } else {
        Vec::new()
    };
    let mut trace_f = if traceback {
        vec![TRACE_F_FROM_H; (n + 1) * (m + 1)]
    } else {
        Vec::new()
    };

    h_row[0] = 0.0;
    if traceback {
        trace_h[0] = DIR_DIAG;
    }
    for j in 1..=m {
        h_row[j] = end_gap_open + end_gap_extend * (j as f32 - 1.0);
        if traceback {
            trace_h[j] = DIR_INS;
        }
    }

    for i in 1..=n {
        let t = target.codes[i - 1];
        let mut h_diag = h_row[0];
        h_row[0] = end_gap_open + end_gap_extend * (i as f32 - 1.0);
        if traceback {
            trace_h[i * (m + 1)] = DIR_DEL;
        }
        let mut f = neg_inf;
        for j in 1..=m {
            let del_gap_o = if scoring.end_gap && j == m {
                end_gap_open
            } else {
                gap_open
            };
            let del_gap_e = if scoring.end_gap && j == m {
                end_gap_extend
            } else {
                gap_extend
            };
            let ins_gap_o = if scoring.end_gap && i == n {
                end_gap_open
            } else {
                gap_open
            };
            let ins_gap_e = if scoring.end_gap && i == n {
                end_gap_extend
            } else {
                gap_extend
            };
            let h_up = h_row[j];
            let e_open = h_up + del_gap_o;
            let e_ext = e_row[j] + del_gap_e;
            let e_from_ext = e_ext > e_open;
            e_row[j] = if e_from_ext { e_ext } else { e_open };
            let f_open = h_row[j - 1] + ins_gap_o;
            let f_ext = f + ins_gap_e;
            let f_from_ext = f_ext > f_open;
            f = if f_from_ext { f_ext } else { f_open };
            let score_diag = h_diag + scoring.score(query.codes[j - 1], t) as f32;
            let mut h = score_diag;
            let mut d = DIR_DIAG;
            // Tie-breaking policy:
            // DIAG > DEL > INS (because we use strict > comparisons)
            // Multiple optimal alignments may exist; this is intentional.
            if e_row[j] > h {
                h = e_row[j];
                d = DIR_DEL;
            }
            if f > h {
                h = f;
                d = DIR_INS;
            }
            if traceback {
                let idx = i * (m + 1) + j;
                // trace_h encodes the predecessor of H(i,j); trace_e/f encode the predecessor of E/F.
                trace_e[idx] = if e_from_ext {
                    TRACE_E_FROM_E
                } else {
                    TRACE_E_FROM_H
                };
                trace_f[idx] = if f_from_ext {
                    TRACE_F_FROM_F
                } else {
                    TRACE_F_FROM_H
                };
                trace_h[idx] = d;
            }
            h_row[j] = h;
            h_diag = h_up;
        }
    }

    let score = h_row[m];
    if !traceback {
        return AlignmentResult {
            score,
            query_end: m - 1,
            target_end: n - 1,
            query_start: Some(0),
            target_start: Some(0),
            cigar: None,
        };
    }

    let mut i = n;
    let mut j = m;
    let mut rev_ops: Vec<(CigarOp, usize)> = Vec::new();
    let mut state = 0u8; // 0=H, 1=E(Del), 2=F(Ins)

    while i > 0 || j > 0 {
        if i == 0 {
            // Boundary: remaining query positions are insertions.
            // These form a single gap block from the row-0 initialization.
            push_rev(&mut rev_ops, CigarOp::Ins, j);
            break;
        }
        if j == 0 {
            push_rev(&mut rev_ops, CigarOp::Del, i);
            break;
        }
        match state {
            0 => {
                let d = trace_h[i * (m + 1) + j];
                match d {
                    DIR_DIAG => {
                        push_rev(&mut rev_ops, CigarOp::Match, 1);
                        i -= 1;
                        j -= 1;
                    }
                    DIR_DEL => {
                        state = 1;
                    }
                    DIR_INS => {
                        state = 2;
                    }
                    _ => break,
                }
            }
            1 => {
                let d = trace_e[i * (m + 1) + j];
                let extending = d == TRACE_E_FROM_E;
                push_rev(&mut rev_ops, CigarOp::Del, 1);
                i -= 1;
                if !extending {
                    state = 0;
                }
            }
            2 => {
                let d = trace_f[i * (m + 1) + j];
                let extending = d == TRACE_F_FROM_F;
                push_rev(&mut rev_ops, CigarOp::Ins, 1);
                j -= 1;
                if !extending {
                    state = 0;
                }
            }
            _ => break,
        }
    }

    AlignmentResult {
        score,
        query_end: m - 1,
        target_end: n - 1,
        query_start: Some(0),
        target_start: Some(0),
        cigar: Some(finalize_cigar(rev_ops)),
    }
}
