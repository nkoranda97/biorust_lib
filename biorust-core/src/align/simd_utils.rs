use super::encode::EncodedSeq;
use super::types::Scoring;
use wide::i16x16;

pub const LANES: usize = 16;

#[inline]
pub fn shift_left(v: i16x16, insert: i16) -> i16x16 {
    let mut arr = v.to_array();
    for i in (1..LANES).rev() {
        arr[i] = arr[i - 1];
    }
    arr[0] = insert;
    i16x16::from(arr)
}

pub fn build_profile(query: &EncodedSeq, scoring: &Scoring) -> Vec<i16x16> {
    let m = query.codes.len();
    let seg_len = m.div_ceil(LANES);
    let alphabet = if scoring.matrix.is_some() {
        scoring
            .alphabet_size
            .expect("alphabet_size must be set when matrix is present")
    } else {
        query.alphabet_size
    };
    let mut profile = vec![i16x16::splat(0); alphabet * seg_len];
    for a in 0..alphabet {
        for seg in 0..seg_len {
            let mut lane_vals = [0i16; LANES];
            for (lane, slot) in lane_vals.iter_mut().enumerate() {
                let idx = lane * seg_len + seg;
                *slot = if idx < m {
                    scoring.score(query.codes[idx], a as u8)
                } else {
                    0
                };
            }
            profile[a * seg_len + seg] = i16x16::from(lane_vals);
        }
    }
    profile
}
