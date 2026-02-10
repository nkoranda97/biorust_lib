use super::encode::EncodedSeq;
use super::simd_utils::{build_profile, shift_left, LANES};
use super::types::Scoring;
use wide::i16x16;

pub fn align_global_score(
    query: &EncodedSeq,
    target: &EncodedSeq,
    scoring: &Scoring,
) -> (f32, usize, usize) {
    let m = query.codes.len();
    let n = target.codes.len();
    if m == 0 || n == 0 {
        let len = if m == 0 { n } else { m };
        let score = if len == 0 {
            0.0
        } else {
            scoring.gap_open + scoring.gap_extend * (len as f32 - 1.0)
        };
        return (score, m.saturating_sub(1), n.saturating_sub(1));
    }

    let seg_len = m.div_ceil(LANES);
    let profile = build_profile(query, scoring);

    let neg_inf = i16::MIN / 2;
    let v_neg_inf = i16x16::splat(neg_inf);
    let gap_open = scoring.gap_open_i16();
    let gap_extend = scoring.gap_extend_i16();
    let v_gap_o = i16x16::splat(gap_open);
    let v_gap_e = i16x16::splat(gap_extend);

    let mut h_prev = vec![v_neg_inf; seg_len];
    let mut h = vec![v_neg_inf; seg_len];
    let mut e = vec![v_neg_inf; seg_len];

    let last_seg = seg_len.saturating_sub(1);
    let needs_mask = m % LANES != 0;
    // Precompute a clamp vector: valid lanes get i16::MAX (no-op under min),
    // invalid lanes get neg_inf. A single v_h.min(v_clamp) replaces the
    // to_array()/from() round-trip per masked segment.
    let v_clamp = if needs_mask {
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

    // Initialize row 0 (i = 0)
    for seg in 0..seg_len {
        let mut lane_vals = [neg_inf; LANES];
        for (lane, slot) in lane_vals.iter_mut().enumerate() {
            let idx = lane * seg_len + seg;
            if idx < m {
                let val = -(gap_open as i32) - (gap_extend as i32) * (idx as i32);
                *slot = val.clamp(neg_inf as i32, i16::MAX as i32) as i16;
            }
        }
        let v = i16x16::from(lane_vals);
        h_prev[seg] = v;
        e[seg] = v - v_gap_o;
    }

    let mut h_left_prev: i16 = 0;

    for (t_idx, &tb) in target.codes.iter().enumerate() {
        let mut v_f = v_neg_inf;
        let mut v_h_diag = shift_left(h_prev[seg_len - 1], h_left_prev);
        let prof_base = tb as usize * seg_len;

        for i in 0..seg_len {
            let v_h_old = h_prev[i];
            let v_p = profile[prof_base + i];
            let v_e = e[i];

            let mut v_h = v_h_diag + v_p;
            v_h = v_h.max(v_e);
            v_h = v_h.max(v_f);

            if i == last_seg && needs_mask {
                v_h = v_h.min(v_clamp);
            }
            h[i] = v_h;

            let v_h_gap = v_h - v_gap_o;
            let mut v_e_new = (v_e - v_gap_e).max(v_h_gap);
            if i == last_seg && needs_mask {
                v_e_new = v_e_new.min(v_clamp);
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
                    v_h_i = v_h_i.min(v_clamp);
                }
                *h_slot = v_h_i;
                let v_h_gap = v_h_i - v_gap_o;
                v_f = (v_f - v_gap_e).max(v_h_gap);
            }
        }

        // update leftmost column for next row
        // Use saturating arithmetic to prevent i16 overflow for long sequences
        let h_left_val = gap_open
            .saturating_neg()
            .saturating_sub(gap_extend.saturating_mul(t_idx as i16));
        h_left_prev = h_left_val;

        std::mem::swap(&mut h_prev, &mut h);
    }

    let last_idx = m - 1;
    let seg = last_idx % seg_len;
    let lane = last_idx / seg_len;
    let score = h_prev[seg].to_array()[lane] as f32;

    (score, m - 1, n - 1)
}
