use crate::error::{BioError, BioResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DnaDistanceModel {
    PDistance,
    JukesCantor,
    Kimura2P,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinDistanceModel {
    PDistance,
    Poisson,
}

#[derive(Debug, Clone)]
pub struct DistanceMatrix {
    labels: Vec<Box<str>>,
    data: Vec<f64>,
    n: usize,
}

impl DistanceMatrix {
    pub fn new(labels: Vec<Box<str>>, data: Vec<f64>) -> Self {
        let n = labels.len();
        assert_eq!(
            data.len(),
            n * n,
            "distance matrix data length mismatch: expected {}, got {}",
            n * n,
            data.len()
        );
        Self { labels, data, n }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn labels(&self) -> &[Box<str>] {
        &self.labels
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    pub fn set(&mut self, i: usize, j: usize, val: f64) {
        self.data[i * self.n + j] = val;
        self.data[j * self.n + i] = val;
    }
}

#[inline]
fn is_gap(b: u8) -> bool {
    b == b'-' || b == b'.'
}

fn count_dna_differences(a: &[u8], b: &[u8]) -> (usize, usize, usize) {
    let mut transitions = 0usize;
    let mut transversions = 0usize;
    let mut valid = 0usize;

    for (&x, &y) in a.iter().zip(b.iter()) {
        if is_gap(x) || is_gap(y) {
            continue;
        }
        let xu = x.to_ascii_uppercase();
        let yu = y.to_ascii_uppercase();
        if !matches!(xu, b'A' | b'C' | b'G' | b'T') || !matches!(yu, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        valid += 1;
        if xu == yu {
            continue;
        }
        // Transitions: A<->G, C<->T
        let is_ts = matches!(
            (xu, yu),
            (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
        );
        if is_ts {
            transitions += 1;
        } else {
            transversions += 1;
        }
    }

    (transitions, transversions, valid)
}

fn count_protein_differences(a: &[u8], b: &[u8]) -> (usize, usize) {
    let mut mismatches = 0usize;
    let mut valid = 0usize;

    for (&x, &y) in a.iter().zip(b.iter()) {
        if is_gap(x) || is_gap(y) {
            continue;
        }
        valid += 1;
        if !x.eq_ignore_ascii_case(&y) {
            mismatches += 1;
        }
    }

    (mismatches, valid)
}

fn compute_dna_pair_distance(
    a: &[u8],
    b: &[u8],
    model: DnaDistanceModel,
    i: usize,
    j: usize,
) -> BioResult<f64> {
    let (ts, tv, valid) = count_dna_differences(a, b);
    if valid == 0 {
        return Err(BioError::NoValidSites { i, j });
    }

    match model {
        DnaDistanceModel::PDistance => Ok((ts + tv) as f64 / valid as f64),
        DnaDistanceModel::JukesCantor => {
            let p = (ts + tv) as f64 / valid as f64;
            let arg = 1.0 - 4.0 * p / 3.0;
            if arg <= 0.0 {
                return Err(BioError::SaturatedDistance {
                    i,
                    j,
                    model: "JukesCantor".into(),
                });
            }
            Ok(-0.75 * arg.ln())
        }
        DnaDistanceModel::Kimura2P => {
            let p = ts as f64 / valid as f64;
            let q = tv as f64 / valid as f64;
            let a1 = 1.0 - 2.0 * p - q;
            let a2 = 1.0 - 2.0 * q;
            if a1 <= 0.0 || a2 <= 0.0 {
                return Err(BioError::SaturatedDistance {
                    i,
                    j,
                    model: "Kimura2P".into(),
                });
            }
            Ok(-0.5 * a1.ln() - 0.25 * a2.ln())
        }
    }
}

fn compute_protein_pair_distance(
    a: &[u8],
    b: &[u8],
    model: ProteinDistanceModel,
    i: usize,
    j: usize,
) -> BioResult<f64> {
    let (mismatches, valid) = count_protein_differences(a, b);
    if valid == 0 {
        return Err(BioError::NoValidSites { i, j });
    }

    let p = mismatches as f64 / valid as f64;

    match model {
        ProteinDistanceModel::PDistance => Ok(p),
        ProteinDistanceModel::Poisson => {
            let arg = 1.0 - p;
            if arg <= 0.0 {
                return Err(BioError::SaturatedDistance {
                    i,
                    j,
                    model: "Poisson".into(),
                });
            }
            Ok(-arg.ln())
        }
    }
}

fn validate_distance_inputs(seqs: &[&[u8]], labels: &[Box<str>]) -> BioResult<()> {
    let n = seqs.len();
    if n < 2 {
        return Err(BioError::TooFewSequences { n });
    }
    if labels.len() != n {
        return Err(BioError::LabelCountMismatch {
            labels: labels.len(),
            seqs: n,
        });
    }
    let expected_len = seqs[0].len();
    for (idx, seq) in seqs.iter().enumerate() {
        if seq.len() != expected_len {
            return Err(BioError::SequenceLengthMismatch {
                index: idx,
                len: seq.len(),
                expected: expected_len,
            });
        }
    }
    Ok(())
}

pub fn dna_distance_matrix(
    seqs: &[&[u8]],
    labels: Vec<Box<str>>,
    model: DnaDistanceModel,
) -> BioResult<DistanceMatrix> {
    validate_distance_inputs(seqs, &labels)?;
    let n = seqs.len();

    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .collect();

    let results: BioResult<Vec<(usize, usize, f64)>> = par_try_map!(&pairs, |&(i, j)| {
        compute_dna_pair_distance(seqs[i], seqs[j], model, i, j).map(|d| (i, j, d))
    });

    let mut data = vec![0.0f64; n * n];
    for (i, j, d) in results? {
        data[i * n + j] = d;
        data[j * n + i] = d;
    }

    Ok(DistanceMatrix::new(labels, data))
}

pub fn protein_distance_matrix(
    seqs: &[&[u8]],
    labels: Vec<Box<str>>,
    model: ProteinDistanceModel,
) -> BioResult<DistanceMatrix> {
    validate_distance_inputs(seqs, &labels)?;
    let n = seqs.len();

    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .collect();

    let results: BioResult<Vec<(usize, usize, f64)>> = par_try_map!(&pairs, |&(i, j)| {
        compute_protein_pair_distance(seqs[i], seqs[j], model, i, j).map(|d| (i, j, d))
    });

    let mut data = vec![0.0f64; n * n];
    for (i, j, d) in results? {
        data[i * n + j] = d;
        data[j * n + i] = d;
    }

    Ok(DistanceMatrix::new(labels, data))
}
