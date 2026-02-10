use super::encode::EncodedSeq;
use super::simd_utils::{build_profile, shift_left, LANES};
use super::types::Scoring;
use wide::i16x16;

pub fn align_local_score(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
) -> (f32, usize, usize) {
    let m = query.codes.len();
    let n = target.codes.len();
    if m == 0 || n == 0 {
        return (0.0, 0, 0);
    }

    let seg_len = m.div_ceil(LANES);
    let profile = build_profile(query, scoring);

    let neg_inf = i16::MIN / 2;
    let v_zero = i16x16::splat(0);
    let v_neg_inf = i16x16::splat(neg_inf);
    // SIMD kernel expects positive gap penalties (subtracted from scores).
    let gap_open = scoring.gap_open_i16();
    let gap_extend = scoring.gap_extend_i16();
    let v_gap_o = i16x16::splat(gap_open);
    let v_gap_e = i16x16::splat(gap_extend);

    let mut h_prev = vec![v_zero; seg_len];
    let mut h = vec![v_zero; seg_len];
    let mut e = vec![v_neg_inf; seg_len];

    let last_seg = seg_len.saturating_sub(1);
    let needs_mask = m % LANES != 0;
    // Precompute clamp vectors for invalid lanes.
    // H invalid lanes → 0 (local alignment floors at zero).
    // E invalid lanes → neg_inf (gap scores stay suppressed).
    let v_clamp_h = if needs_mask {
        let mut arr = [i16::MAX; LANES];
        for (lane, slot) in arr.iter_mut().enumerate() {
            let idx = lane * seg_len + last_seg;
            if idx >= m {
                *slot = 0;
            }
        }
        i16x16::from(arr)
    } else {
        i16x16::splat(i16::MAX)
    };
    let v_clamp_e = if needs_mask {
        let mut arr = [i16::MAX; LANES];
        for (lane, slot) in arr.iter_mut().enumerate() {
            let idx = lane * seg_len + last_seg;
            if idx >= m {
                *slot = neg_inf;
            }
        }
        i16x16::from(arr)
    } else {
        i16x16::splat(i16::MAX)
    };

    let mut max_score: i16 = 0;
    let mut end_q: usize = 0;
    let mut end_t: usize = 0;

    for (t_idx, &tb) in target.codes.iter().enumerate() {
        let mut v_f = v_neg_inf;
        let mut v_h_diag = shift_left(h_prev[seg_len - 1], 0);
        let prof_base = tb as usize * seg_len;

        for i in 0..seg_len {
            let v_h_old = h_prev[i];
            let v_p = profile[prof_base + i];
            let v_e = e[i];

            let mut v_h = v_h_diag + v_p;
            v_h = v_h.max(v_e);
            v_h = v_h.max(v_f);
            v_h = v_h.max(v_zero);

            if i == last_seg && needs_mask {
                v_h = v_h.min(v_clamp_h);
            }
            h[i] = v_h;

            let v_h_gap = v_h - v_gap_o;
            let mut v_e_new = (v_e - v_gap_e).max(v_h_gap);
            if i == last_seg && needs_mask {
                v_e_new = v_e_new.min(v_clamp_e);
            }
            e[i] = v_e_new;
            v_f = (v_f - v_gap_e).max(v_h_gap);

            v_h_diag = v_h_old;
        }

        // Lazy F loop
        for _ in 0..LANES {
            v_f = shift_left(v_f, neg_inf);
            for (i, h_slot) in h.iter_mut().enumerate() {
                let mut v_h_i = (*h_slot).max(v_f);
                if i == last_seg && needs_mask {
                    v_h_i = v_h_i.min(v_clamp_h);
                }
                *h_slot = v_h_i;
                let v_h_gap = v_h_i - v_gap_o;
                v_f = (v_f - v_gap_e).max(v_h_gap);
            }
            let any_pos = v_f.to_array().iter().take(LANES).any(|&v| v > 0);
            if !any_pos {
                break;
            }
        }

        // Track max score and end position
        for (i, h_vec) in h.iter().enumerate() {
            let arr = h_vec.to_array();
            for (lane, &val) in arr.iter().enumerate().take(LANES) {
                let q_idx = lane * seg_len + i;
                if q_idx >= m {
                    continue;
                }
                if val > max_score {
                    max_score = val;
                    end_q = q_idx;
                    end_t = t_idx;
                }
            }
        }

        std::mem::swap(&mut h_prev, &mut h);
    }

    (max_score as f32, end_q, end_t)
}

#[cfg(test)]
pub(crate) fn align_local_score_rows(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
) -> Vec<i32> {
    let m = query.codes.len();
    let n = target.codes.len();
    if m == 0 || n == 0 {
        return vec![0; n];
    }

    let seg_len = m.div_ceil(LANES);
    let profile = build_profile(query, scoring);

    let neg_inf = i16::MIN / 2;
    let v_zero = i16x16::splat(0);
    let v_neg_inf = i16x16::splat(neg_inf);
    let gap_open = scoring.gap_open_i16();
    let gap_extend = scoring.gap_extend_i16();
    let v_gap_o = i16x16::splat(gap_open);
    let v_gap_e = i16x16::splat(gap_extend);

    let mut h_prev = vec![v_zero; seg_len];
    let mut h = vec![v_zero; seg_len];
    let mut e = vec![v_neg_inf; seg_len];

    let last_seg = seg_len.saturating_sub(1);
    let needs_mask = m % LANES != 0;
    let v_clamp_h = if needs_mask {
        let mut arr = [i16::MAX; LANES];
        for (lane, slot) in arr.iter_mut().enumerate() {
            let idx = lane * seg_len + last_seg;
            if idx >= m {
                *slot = 0;
            }
        }
        i16x16::from(arr)
    } else {
        i16x16::splat(i16::MAX)
    };
    let v_clamp_e = if needs_mask {
        let mut arr = [i16::MAX; LANES];
        for (lane, slot) in arr.iter_mut().enumerate() {
            let idx = lane * seg_len + last_seg;
            if idx >= m {
                *slot = neg_inf;
            }
        }
        i16x16::from(arr)
    } else {
        i16x16::splat(i16::MAX)
    };

    let mut rows = Vec::with_capacity(n);

    for &tb in &target.codes {
        let mut v_f = v_neg_inf;
        let mut v_h_diag = shift_left(h_prev[seg_len - 1], 0);
        let prof_base = tb as usize * seg_len;

        for i in 0..seg_len {
            let v_h_old = h_prev[i];
            let v_p = profile[prof_base + i];
            let v_e = e[i];

            let mut v_h = v_h_diag + v_p;
            v_h = v_h.max(v_e);
            v_h = v_h.max(v_f);
            v_h = v_h.max(v_zero);

            if i == last_seg && needs_mask {
                v_h = v_h.min(v_clamp_h);
            }
            h[i] = v_h;

            let v_h_gap = v_h - v_gap_o;
            let mut v_e_new = (v_e - v_gap_e).max(v_h_gap);
            if i == last_seg && needs_mask {
                v_e_new = v_e_new.min(v_clamp_e);
            }
            e[i] = v_e_new;
            v_f = (v_f - v_gap_e).max(v_h_gap);

            v_h_diag = v_h_old;
        }

        for _ in 0..LANES {
            v_f = shift_left(v_f, neg_inf);
            for (i, h_slot) in h.iter_mut().enumerate() {
                let mut v_h_i = (*h_slot).max(v_f);
                if i == last_seg && needs_mask {
                    v_h_i = v_h_i.min(v_clamp_h);
                }
                *h_slot = v_h_i;
                let v_h_gap = v_h_i - v_gap_o;
                v_f = (v_f - v_gap_e).max(v_h_gap);
            }
            let any_pos = v_f.to_array().iter().take(LANES).any(|&v| v > 0);
            if !any_pos {
                break;
            }
        }

        let mut row_max = 0i16;
        for (i, h_vec) in h.iter().enumerate() {
            let arr = h_vec.to_array();
            for (lane, &val) in arr.iter().enumerate().take(LANES) {
                let q_idx = lane * seg_len + i;
                if q_idx >= m {
                    continue;
                }
                if val > row_max {
                    row_max = val;
                }
            }
        }
        rows.push(row_max as i32);

        std::mem::swap(&mut h_prev, &mut h);
    }

    rows
}
