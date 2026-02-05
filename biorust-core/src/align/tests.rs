use super::encode::encode_dna;
#[cfg(feature = "simd")]
use super::local_simd::align_local_score;
use super::scalar_ref::{align_global_scalar, align_local_scalar};
use super::types::{Cigar, CigarOp, Scoring};
use super::{align_global, align_local, score_alignment_from_cigar};

#[cfg(feature = "simd")]
use proptest::prelude::*;

const GAP: u8 = 0xFF;

fn cigar_consumed_lengths(cigar: &Cigar) -> (usize, usize) {
    let mut q_len = 0usize;
    let mut t_len = 0usize;
    for (op, len) in &cigar.ops {
        match op {
            CigarOp::Match => {
                q_len += *len;
                t_len += *len;
            }
            CigarOp::Ins => {
                q_len += *len;
            }
            CigarOp::Del => {
                t_len += *len;
            }
        }
    }
    (q_len, t_len)
}

fn aligned_from_cigar(query: &[u8], target: &[u8], cigar: &Cigar) -> (Vec<u8>, Vec<u8>) {
    let mut aligned_q = Vec::new();
    let mut aligned_t = Vec::new();
    let mut qi = 0usize;
    let mut ti = 0usize;
    for (op, len) in &cigar.ops {
        match op {
            CigarOp::Match => {
                for _ in 0..*len {
                    aligned_q.push(query[qi]);
                    aligned_t.push(target[ti]);
                    qi += 1;
                    ti += 1;
                }
            }
            CigarOp::Ins => {
                for _ in 0..*len {
                    aligned_q.push(query[qi]);
                    aligned_t.push(GAP);
                    qi += 1;
                }
            }
            CigarOp::Del => {
                for _ in 0..*len {
                    aligned_q.push(GAP);
                    aligned_t.push(target[ti]);
                    ti += 1;
                }
            }
        }
    }
    (aligned_q, aligned_t)
}

fn assert_ungap_roundtrip(
    query: &[u8],
    target: &[u8],
    cigar: &Cigar,
    q_start: usize,
    t_start: usize,
) {
    let (q_len, t_len) = cigar_consumed_lengths(cigar);
    let q_slice = &query[q_start..q_start + q_len];
    let t_slice = &target[t_start..t_start + t_len];
    let (aligned_q, aligned_t) = aligned_from_cigar(q_slice, t_slice, cigar);
    let ungapped_q: Vec<u8> = aligned_q.iter().copied().filter(|&b| b != GAP).collect();
    let ungapped_t: Vec<u8> = aligned_t.iter().copied().filter(|&b| b != GAP).collect();
    assert_eq!(ungapped_q, q_slice);
    assert_eq!(ungapped_t, t_slice);
    assert_eq!(aligned_q.len(), aligned_t.len());
}

fn rescore_from_cigar(
    query: &[u8],
    target: &[u8],
    cigar: &Cigar,
    scoring: &Scoring,
    q_start: usize,
    t_start: usize,
) -> f32 {
    let q_enc = encode_dna(query).unwrap();
    let t_enc = encode_dna(target).unwrap();
    let (q_len, t_len) = cigar_consumed_lengths(cigar);
    let q_slice = &q_enc.codes[q_start..q_start + q_len];
    let t_slice = &t_enc.codes[t_start..t_start + t_len];
    score_alignment_from_cigar(q_slice, t_slice, cigar, scoring)
}

#[test]
fn encode_dna_valid() {
    let enc = encode_dna(b"ACGTNacgtn").unwrap();
    assert_eq!(enc.alphabet_size, 15);
    assert_eq!(enc.codes.len(), 10);
}

#[test]
fn local_scalar_simple_match() {
    let q = encode_dna(b"ACGT").unwrap();
    let t = encode_dna(b"ACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_local_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, 8.0);
    assert_eq!(res.query_end, 3);
    assert_eq!(res.target_end, 3);
    assert_eq!(res.query_start, Some(0));
    assert_eq!(res.target_start, Some(0));
    assert!(res.cigar.is_some());
}

#[test]
fn global_scalar_simple_match() {
    let q = encode_dna(b"ACGT").unwrap();
    let t = encode_dna(b"ACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, 8.0);
    assert_eq!(res.query_end, 3);
    assert_eq!(res.target_end, 3);
    assert!(res.cigar.is_some());
}

#[test]
fn global_scalar_gap() {
    let q = encode_dna(b"ACGT").unwrap();
    let t = encode_dna(b"ACG").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, 4.0);
}

#[test]
fn align_local_with_traceback() {
    let q = encode_dna(b"ACGTAC").unwrap();
    let t = encode_dna(b"TTACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_local(&q, &t, &scoring, true);
    assert!(res.cigar.is_some());
    assert!(res.query_start.is_some());
    assert!(res.target_start.is_some());
}

#[test]
fn align_global_with_traceback() {
    let q = encode_dna(b"ACGT").unwrap();
    let t = encode_dna(b"ACG").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global(&q, &t, &scoring, true);
    assert!(res.cigar.is_some());
    assert_eq!(res.query_start, Some(0));
    assert_eq!(res.target_start, Some(0));
}

#[test]
fn global_ungap_roundtrip() {
    let q_bytes = b"ACGTAC";
    let t_bytes = b"ACG";
    let q = encode_dna(q_bytes).unwrap();
    let t = encode_dna(t_bytes).unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    let cigar = res.cigar.as_ref().unwrap();
    assert_ungap_roundtrip(&q.codes, &t.codes, cigar, 0, 0);
}

#[test]
fn global_rescore_matches_dp() {
    let q_bytes = b"ACGTACGT";
    let t_bytes = b"ACGT";
    let q = encode_dna(q_bytes).unwrap();
    let t = encode_dna(t_bytes).unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    let cigar = res.cigar.as_ref().unwrap();
    let recomputed = rescore_from_cigar(q_bytes, t_bytes, cigar, &scoring, 0, 0);
    assert_eq!(res.score, recomputed);
}

#[test]
fn local_rescore_matches_dp() {
    let q_bytes = b"NNACGTNN";
    let t_bytes = b"ACGT";
    let q = encode_dna(q_bytes).unwrap();
    let t = encode_dna(t_bytes).unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_local_scalar(&q, &t, &scoring, true);
    let cigar = res.cigar.as_ref().unwrap();
    let qs = res.query_start.unwrap();
    let ts = res.target_start.unwrap();
    let recomputed = rescore_from_cigar(q_bytes, t_bytes, cigar, &scoring, qs, ts);
    assert_eq!(res.score, recomputed);
}

#[cfg(feature = "simd")]
proptest! {
    #[test]
    fn simd_scalar_parity_global(
        q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..40),
        t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..40),
    ) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(2, -1, -2.0, -1.0);
        prop_assume!(scoring.simd_compatible());
        let scalar = align_global_scalar(&q_enc, &t_enc, &scoring, false);
        let (simd_score, _, _) = super::global_simd::align_global_score(&q_enc, &t_enc, &scoring);
        prop_assert_eq!(simd_score, scalar.score);
    }
}

#[cfg(feature = "simd")]
proptest! {
    #[test]
    fn local_simd_matches_scalar(q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..40),
                                  t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..40)) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(2, -1, -2.0, -1.0);
        let scalar = align_local_scalar(&q_enc, &t_enc, &scoring, false);
        let (simd_score, _, _) = align_local_score(&q_enc, &t_enc, &scoring);
        prop_assert_eq!(simd_score, scalar.score);
    }
}

#[cfg(feature = "simd")]
proptest! {
    #[test]
    fn global_simd_matches_scalar(q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..30),
                                  t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..30)) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(2, -1, -2.0, -1.0);
        let scalar = align_global_scalar(&q_enc, &t_enc, &scoring, false);
        let (simd_score, _, _) = super::global_simd::align_global_score(&q_enc, &t_enc, &scoring);
        prop_assert_eq!(simd_score, scalar.score);
    }
}

#[cfg(feature = "simd")]
#[test]
fn local_simd_debug_rows() {
    let q = b"CGTTACCGGAAG";
    let t = b"AGATCCCAAG";
    let q_enc = encode_dna(q).unwrap();
    let t_enc = encode_dna(t).unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let simd_rows = super::local_simd::align_local_score_rows(&q_enc, &t_enc, &scoring);
    let mut scalar_rows = Vec::new();
    let m = q_enc.codes.len();
    let n = t_enc.codes.len();
    let neg_inf: i32 = i32::MIN / 4;
    let gap_open = scoring.gap_open.round() as i32;
    let gap_extend = scoring.gap_extend.round() as i32;
    let mut h_row = vec![0i32; m + 1];
    let mut e_row = vec![neg_inf; m + 1];
    for i in 1..=n {
        let tbase = t_enc.codes[i - 1];
        let mut h_diag = 0i32;
        let mut f = neg_inf;
        h_row[0] = 0;
        for j in 1..=m {
            let h_up = h_row[j];
            e_row[j] = (h_up + gap_open).max(e_row[j] + gap_extend);
            f = (h_row[j - 1] + gap_open).max(f + gap_extend);
            let score_diag = h_diag + scoring.score(q_enc.codes[j - 1], tbase) as i32;
            let mut h = score_diag.max(e_row[j]).max(f);
            if h < 0 {
                h = 0;
            }
            h_row[j] = h;
            h_diag = h_up;
        }
        let row_max = *h_row.iter().max().unwrap();
        scalar_rows.push(row_max);
    }
    if simd_rows != scalar_rows {
        eprintln!("simd_rows:   {:?}", simd_rows);
        eprintln!("scalar_rows: {:?}", scalar_rows);
    }
    assert_eq!(simd_rows, scalar_rows);
}

// ---- Known-answer tests ----

#[test]
fn global_known_answer_insertion() {
    // Query: ACGT, Target: AGT => one deletion from target perspective (or insertion in query)
    let q = encode_dna(b"ACGT").unwrap();
    let t = encode_dna(b"AGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    let cigar = res.cigar.as_ref().unwrap();
    let recomputed = rescore_from_cigar(b"ACGT", b"AGT", cigar, &scoring, 0, 0);
    assert_eq!(res.score, recomputed, "cigar: {:?}", cigar.ops);
}

#[test]
fn global_known_answer_deletion() {
    // Query: AGT, Target: ACGT
    let q = encode_dna(b"AGT").unwrap();
    let t = encode_dna(b"ACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    let cigar = res.cigar.as_ref().unwrap();
    let recomputed = rescore_from_cigar(b"AGT", b"ACGT", cigar, &scoring, 0, 0);
    assert_eq!(res.score, recomputed, "cigar: {:?}", cigar.ops);
}

#[test]
fn global_known_answer_all_mismatches() {
    let q = encode_dna(b"AAAA").unwrap();
    let t = encode_dna(b"CCCC").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    // All mismatches: 4 * -1 = -4
    assert_eq!(res.score, -4.0);
    let cigar = res.cigar.as_ref().unwrap();
    let recomputed = rescore_from_cigar(b"AAAA", b"CCCC", cigar, &scoring, 0, 0);
    assert_eq!(res.score, recomputed);
}

#[test]
fn global_known_answer_long_gap() {
    // Query has extra bases that must be gapped
    let q = encode_dna(b"ACGTACGT").unwrap();
    let t = encode_dna(b"ACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    let cigar = res.cigar.as_ref().unwrap();
    let recomputed = rescore_from_cigar(b"ACGTACGT", b"ACGT", cigar, &scoring, 0, 0);
    assert_eq!(res.score, recomputed, "cigar: {:?}", cigar.ops);
}

#[test]
fn local_known_answer_partial_match() {
    let q = encode_dna(b"NNACGTNN").unwrap();
    let t = encode_dna(b"ACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_local_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, 8.0); // perfect 4-match
    let cigar = res.cigar.as_ref().unwrap();
    let qs = res.query_start.unwrap();
    let ts = res.target_start.unwrap();
    let recomputed = rescore_from_cigar(b"NNACGTNN", b"ACGT", cigar, &scoring, qs, ts);
    assert_eq!(res.score, recomputed, "cigar: {:?}", cigar.ops);
}

#[test]
fn global_single_char_sequences() {
    let q = encode_dna(b"A").unwrap();
    let t = encode_dna(b"A").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, 2.0);
    let cigar = res.cigar.as_ref().unwrap();
    assert_eq!(cigar.ops, vec![(CigarOp::Match, 1)]);
}

#[test]
fn global_single_char_mismatch() {
    let q = encode_dna(b"A").unwrap();
    let t = encode_dna(b"C").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, -1.0);
    let cigar = res.cigar.as_ref().unwrap();
    assert_eq!(cigar.ops, vec![(CigarOp::Match, 1)]);
}

#[test]
fn global_all_gaps_query_empty() {
    let q = encode_dna(b"").unwrap();
    let t = encode_dna(b"ACGT").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    // gap_open + gap_extend * 3 = -2 + (-1)*3 = -5
    assert_eq!(res.score, -5.0);
}

#[test]
fn global_all_gaps_target_empty() {
    let q = encode_dna(b"ACGT").unwrap();
    let t = encode_dna(b"").unwrap();
    let scoring = Scoring::simple(2, -1, -2.0, -1.0);
    let res = align_global_scalar(&q, &t, &scoring, true);
    assert_eq!(res.score, -5.0);
}

// ---- Traceback validation with varied scoring ----

#[test]
fn global_traceback_varied_scoring() {
    let cases: Vec<(&[u8], &[u8], Scoring)> = vec![
        (b"ACGTACGT", b"ACGTACGT", Scoring::simple(1, -3, -5.0, -2.0)),
        (b"AAAAAA", b"AAA", Scoring::simple(1, -1, -3.0, -1.0)),
        (b"ACGT", b"TGCA", Scoring::simple(2, -2, -4.0, -1.0)),
        (b"ACGTACGT", b"ACTACG", Scoring::simple(3, -1, -5.0, -2.0)),
        (b"A", b"ACGT", Scoring::simple(1, -1, -2.0, -1.0)),
    ];
    for (q_bytes, t_bytes, scoring) in &cases {
        let q = encode_dna(q_bytes).unwrap();
        let t = encode_dna(t_bytes).unwrap();
        let res = align_global_scalar(&q, &t, scoring, true);
        let cigar = res.cigar.as_ref().unwrap();
        let recomputed = rescore_from_cigar(q_bytes, t_bytes, cigar, scoring, 0, 0);
        assert_eq!(
            res.score,
            recomputed,
            "global q={} t={} cigar={:?}",
            std::str::from_utf8(q_bytes).unwrap(),
            std::str::from_utf8(t_bytes).unwrap(),
            cigar.ops
        );
    }
}

#[test]
fn local_traceback_varied_scoring() {
    let cases: Vec<(&[u8], &[u8], Scoring)> = vec![
        (b"ACGTACGT", b"ACGTACGT", Scoring::simple(1, -3, -5.0, -2.0)),
        (b"AAACCCTTT", b"AAATTT", Scoring::simple(1, -1, -3.0, -1.0)),
        (b"ACGT", b"TGCA", Scoring::simple(2, -2, -4.0, -1.0)),
        (b"GATTACA", b"GCATGCA", Scoring::simple(3, -1, -5.0, -2.0)),
    ];
    for (q_bytes, t_bytes, scoring) in &cases {
        let q = encode_dna(q_bytes).unwrap();
        let t = encode_dna(t_bytes).unwrap();
        let res = align_local_scalar(&q, &t, scoring, true);
        if res.score > 0.0 {
            let cigar = res.cigar.as_ref().unwrap();
            let qs = res.query_start.unwrap();
            let ts = res.target_start.unwrap();
            let recomputed = rescore_from_cigar(q_bytes, t_bytes, cigar, scoring, qs, ts);
            assert_eq!(
                res.score,
                recomputed,
                "local q={} t={} cigar={:?} qs={} ts={}",
                std::str::from_utf8(q_bytes).unwrap(),
                std::str::from_utf8(t_bytes).unwrap(),
                cigar.ops,
                qs,
                ts
            );
        }
    }
}

// ---- Property tests with varied scoring ----

#[cfg(feature = "simd")]
proptest! {
    #[test]
    fn global_simd_matches_scalar_varied(
        q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..25),
        t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..25),
        match_score in 1i16..=4,
        mismatch in -4i16..=-1,
        gap_open in -6i16..=-2,
        gap_ext in -3i16..=-1,
    ) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(match_score, mismatch, gap_open as f32, gap_ext as f32);
        let scalar = align_global_scalar(&q_enc, &t_enc, &scoring, false);
        let (simd_score, _, _) = super::global_simd::align_global_score(&q_enc, &t_enc, &scoring);
        prop_assert_eq!(simd_score, scalar.score,
            "q={:?} t={:?} scoring=({},{},{},{})",
            std::str::from_utf8(&q).unwrap(),
            std::str::from_utf8(&t).unwrap(),
            match_score, mismatch, gap_open, gap_ext);
    }
}

#[cfg(feature = "simd")]
proptest! {
    #[test]
    fn local_simd_matches_scalar_varied(
        q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..25),
        t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..25),
        match_score in 1i16..=4,
        mismatch in -4i16..=-1,
        gap_open in -6i16..=-2,
        gap_ext in -3i16..=-1,
    ) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(match_score, mismatch, gap_open as f32, gap_ext as f32);
        let scalar = align_local_scalar(&q_enc, &t_enc, &scoring, false);
        let (simd_score, _, _) = align_local_score(&q_enc, &t_enc, &scoring);
        prop_assert_eq!(simd_score, scalar.score,
            "q={:?} t={:?} scoring=({},{},{},{})",
            std::str::from_utf8(&q).unwrap(),
            std::str::from_utf8(&t).unwrap(),
            match_score, mismatch, gap_open, gap_ext);
    }
}

// ---- Traceback score recomputation property tests ----

proptest! {
    #[test]
    fn global_traceback_score_matches(
        q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..20),
        t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..20),
    ) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(2, -1, -2.0, -1.0);
        let res = align_global_scalar(&q_enc, &t_enc, &scoring, true);
        let cigar = res.cigar.as_ref().unwrap();
        let recomputed = rescore_from_cigar(&q, &t, cigar, &scoring, 0, 0);
        prop_assert_eq!(res.score, recomputed,
            "q={:?} t={:?} cigar={:?}",
            std::str::from_utf8(&q).unwrap(),
            std::str::from_utf8(&t).unwrap(),
            cigar.ops);
        assert_ungap_roundtrip(&q_enc.codes, &t_enc.codes, cigar, 0, 0);
    }
}

proptest! {
    #[test]
    fn local_traceback_score_matches(
        q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..20),
        t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..20),
    ) {
        let q_enc = encode_dna(&q).unwrap();
        let t_enc = encode_dna(&t).unwrap();
        let scoring = Scoring::simple(2, -1, -2.0, -1.0);
        let res = align_local_scalar(&q_enc, &t_enc, &scoring, true);
        if res.score > 0.0 {
            let cigar = res.cigar.as_ref().unwrap();
            let qs = res.query_start.unwrap();
            let ts = res.target_start.unwrap();
            let recomputed = rescore_from_cigar(&q, &t, cigar, &scoring, qs, ts);
            prop_assert_eq!(res.score, recomputed,
                "q={:?} t={:?} cigar={:?} qs={} ts={}",
                std::str::from_utf8(&q).unwrap(),
                std::str::from_utf8(&t).unwrap(),
                cigar.ops, qs, ts);
            assert_ungap_roundtrip(&q_enc.codes, &t_enc.codes, cigar, qs, ts);
        }
    }
}

proptest! {

    #[test]

    fn global_traceback_score_matches_varied(

        q in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..15),

        t in prop::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..15),

        match_score in 1i16..=4,

        mismatch in -4i16..=-1,

        gap_open in -6i16..=-2,

        gap_ext in -3i16..=-1,

    ) {

        prop_assume!(gap_ext >= gap_open);

        let q_enc = encode_dna(&q).unwrap();

        let t_enc = encode_dna(&t).unwrap();

        let scoring = Scoring::simple(match_score, mismatch, gap_open as f32, gap_ext as f32);

        let res = align_global_scalar(&q_enc, &t_enc, &scoring, true);

        let cigar = res.cigar.as_ref().unwrap();

        let recomputed = rescore_from_cigar(&q, &t, cigar, &scoring, 0, 0);

        prop_assert_eq!(res.score, recomputed,

            "q={:?} t={:?} scoring=({},{},{},{}) cigar={:?}",

            std::str::from_utf8(&q).unwrap(),

            std::str::from_utf8(&t).unwrap(),

            match_score, mismatch, gap_open, gap_ext,

            cigar.ops);

    }

}
