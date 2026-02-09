use crate::alphabets::protein;
use crate::error::{BioError, BioResult};
use crate::seq::bytes::{self, IntoNeedle, Needle};
use crate::seq::traits::SeqBytes;
use std::sync::LazyLock;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ProteinSeq {
    bytes: Vec<u8>,
}

impl ProteinSeq {
    pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
        let alphabet = protein::iupac_alphabet();
        for (pos, &b) in bytes.iter().enumerate() {
            if !alphabet.symbols.contains(b as usize) {
                return Err(BioError::InvalidChar { ch: b as char, pos });
            }
        }
        Ok(Self { bytes })
    }

    #[inline]
    pub(crate) fn from_bytes_unchecked(bytes: Vec<u8>) -> Self {
        Self { bytes }
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    pub fn to_string(&self) -> BioResult<String> {
        // All valid IUPAC protein bytes are valid UTF-8
        std::str::from_utf8(self.as_bytes())
            .map(|s| s.to_string())
            .map_err(|_| BioError::InvalidChar {
                ch: '\u{FFFD}',
                pos: 0,
            })
    }

    pub fn reverse(&self) -> Self {
        let mut out = self.bytes.clone();
        out.reverse();
        Self { bytes: out }
    }

    pub fn count<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::count(self.as_bytes(), needle))
    }

    pub fn count_overlap<'a, N>(&'a self, sub: N) -> BioResult<usize>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::count_overlap(self.as_bytes(), needle))
    }

    pub fn contains<'a, N>(&'a self, sub: N) -> BioResult<bool>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::contains(self.as_bytes(), needle))
    }

    pub fn find<'a, N>(&'a self, sub: N, start: usize, end: usize) -> BioResult<Option<usize>>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::find(self.as_bytes(), needle, start, end))
    }

    pub fn rfind<'a, N>(&'a self, sub: N, start: usize, end: usize) -> BioResult<Option<usize>>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::rfind(self.as_bytes(), needle, start, end))
    }

    pub fn starts_with<'a, N>(&'a self, sub: N) -> BioResult<bool>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::find(self.as_bytes(), needle, 0, self.len()) == Some(0))
    }

    pub fn ends_with<'a, N>(&'a self, sub: N) -> BioResult<bool>
    where
        N: IntoNeedle<'a>,
    {
        let needle = sub.into_needle()?;
        Ok(bytes::rfind(self.as_bytes(), needle, 0, self.len())
            .map(|pos| pos + needle_len(needle) == self.len())
            .unwrap_or(false))
    }

    pub fn counts(&self) -> [u32; 256] {
        let mut counts = [0u32; 256];
        for &b in self.as_bytes() {
            counts[b as usize] += 1;
        }
        counts
    }

    pub fn frequencies(&self) -> [f64; 256] {
        let mut freq = [0.0f64; 256];
        let len = self.len();
        if len == 0 {
            return freq;
        }
        let counts = self.counts();
        let denom = len as f64;
        for (i, &c) in counts.iter().enumerate() {
            if c > 0 {
                freq[i] = c as f64 / denom;
            }
        }
        freq
    }

    pub fn aa_counts_20(&self) -> [u32; 20] {
        let mut counts = [0u32; 20];
        for &b in self.as_bytes() {
            let idx = AA20_INDEX[b as usize];
            if idx >= 0 {
                counts[idx as usize] += 1;
            }
        }
        counts
    }

    pub fn aa_frequencies_20(&self) -> [f64; 20] {
        let mut freq = [0.0f64; 20];
        let len = self.len();
        if len == 0 {
            return freq;
        }
        let counts = self.aa_counts_20();
        let denom = len as f64;
        for (i, &c) in counts.iter().enumerate() {
            if c > 0 {
                freq[i] = c as f64 / denom;
            }
        }
        freq
    }

    pub fn shannon_entropy(&self) -> f64 {
        let len = self.len();
        if len == 0 {
            return 0.0;
        }
        let counts = self.counts();
        let denom = len as f64;
        let mut entropy = 0.0f64;
        for &count in &counts {
            if count == 0 {
                continue;
            }
            let p = count as f64 / denom;
            entropy -= p * p.log2();
        }
        entropy
    }

    pub fn molecular_weight(&self) -> BioResult<f64> {
        if self.is_empty() {
            return Ok(0.0);
        }
        let mut total = 0.0f64;
        for (pos, &b) in self.as_bytes().iter().enumerate() {
            let idx = AA20_INDEX[b as usize];
            if idx < 0 {
                return Err(BioError::InvalidChar { ch: b as char, pos });
            }
            total += AA20_MASS_AVG[idx as usize];
        }
        Ok(total + WATER_MASS)
    }

    pub fn hydrophobicity(&self) -> BioResult<f64> {
        if self.is_empty() {
            return Ok(0.0);
        }
        let mut total = 0.0f64;
        for (pos, &b) in self.as_bytes().iter().enumerate() {
            let idx = AA20_INDEX[b as usize];
            if idx < 0 {
                return Err(BioError::InvalidChar { ch: b as char, pos });
            }
            total += AA20_HYDRO_KD[idx as usize];
        }
        Ok(total / self.len() as f64)
    }

    pub fn hydrophobicity_profile(&self, window: usize) -> BioResult<Vec<f64>> {
        if window == 0 {
            return Err(BioError::InvalidWindow { window });
        }
        let bytes = self.as_bytes();
        if bytes.len() < window {
            return Ok(Vec::new());
        }

        // Pre-map bytes to hydrophobicity values, validating all at once
        let hydro: Vec<f64> = bytes
            .iter()
            .enumerate()
            .map(|(pos, &b)| {
                let idx = AA20_INDEX[b as usize];
                if idx < 0 {
                    Err(BioError::InvalidChar { ch: b as char, pos })
                } else {
                    Ok(AA20_HYDRO_KD[idx as usize])
                }
            })
            .collect::<BioResult<_>>()?;

        // Sliding window: O(n) instead of O(n * window)
        let inv_w = 1.0 / window as f64;
        let mut sum: f64 = hydro[..window].iter().sum();
        let mut out = Vec::with_capacity(bytes.len() - window + 1);
        out.push(sum * inv_w);
        for i in 1..=bytes.len() - window {
            sum += hydro[i + window - 1] - hydro[i - 1];
            out.push(sum * inv_w);
        }
        Ok(out)
    }

    pub fn net_charge(&self, ph: f64) -> BioResult<f64> {
        if self.has_ambiguous() {
            for (pos, &b) in self.as_bytes().iter().enumerate() {
                if AA20_INDEX[b as usize] < 0 {
                    return Err(BioError::InvalidChar { ch: b as char, pos });
                }
            }
        }

        let counts = self.aa_counts_20();
        let n_term = basic_charge(ph, PKA_NTERM);
        let c_term = acidic_charge(ph, PKA_CTERM);
        let mut total = n_term + c_term;

        total += counts[idx('R')] as f64 * basic_charge(ph, PKA_R);
        total += counts[idx('K')] as f64 * basic_charge(ph, PKA_K);
        total += counts[idx('H')] as f64 * basic_charge(ph, PKA_H);
        total += counts[idx('D')] as f64 * acidic_charge(ph, PKA_D);
        total += counts[idx('E')] as f64 * acidic_charge(ph, PKA_E);
        total += counts[idx('C')] as f64 * acidic_charge(ph, PKA_C);
        total += counts[idx('Y')] as f64 * acidic_charge(ph, PKA_Y);

        Ok(total)
    }

    pub fn isoelectric_point(&self) -> BioResult<f64> {
        let mut low = 0.0f64;
        let mut high = 14.0f64;
        for _ in 0..60 {
            let mid = (low + high) / 2.0;
            let charge = self.net_charge(mid)?;
            if charge > 0.0 {
                low = mid;
            } else {
                high = mid;
            }
        }
        Ok((low + high) / 2.0)
    }

    pub fn validate_strict_20(&self) -> BioResult<()> {
        for (pos, &b) in self.as_bytes().iter().enumerate() {
            if AA20_INDEX[b as usize] < 0 {
                return Err(BioError::InvalidChar { ch: b as char, pos });
            }
        }
        Ok(())
    }

    pub fn has_ambiguous(&self) -> bool {
        self.as_bytes().iter().any(|&b| AA20_INDEX[b as usize] < 0)
    }

    pub fn unknown_positions(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for (i, &b) in self.as_bytes().iter().enumerate() {
            if AA20_INDEX[b as usize] < 0 {
                out.push(i);
            }
        }
        out
    }
}

impl SeqBytes for ProteinSeq {
    fn as_bytes(&self) -> &[u8] {
        ProteinSeq::as_bytes(self)
    }

    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self> {
        ProteinSeq::new(bytes)
    }
}

impl<'a> IntoNeedle<'a> for &'a ProteinSeq {
    #[inline]
    fn into_needle(self) -> BioResult<Needle<'a>> {
        Ok(Needle::Bytes(self.as_bytes()))
    }
}

#[inline]
fn needle_len(needle: Needle<'_>) -> usize {
    match needle {
        Needle::Byte(_) => 1,
        Needle::Bytes(bytes) => bytes.len(),
    }
}

static AA20_INDEX: LazyLock<[i8; 256]> = LazyLock::new(|| {
    let mut map = [-1i8; 256];
    for (idx, &b) in AA20.iter().enumerate() {
        map[b as usize] = idx as i8;
        let lower = b.to_ascii_lowercase();
        map[lower as usize] = idx as i8;
    }
    map
});

const AA20: [u8; 20] = *b"ARNDCEQGHILKMFPSTWYV";

const AA20_MASS_AVG: [f64; 20] = [
    71.0788,  // A
    156.1875, // R
    114.1038, // N
    115.0886, // D
    103.1388, // C
    129.1155, // E
    128.1307, // Q
    57.0519,  // G
    137.1411, // H
    113.1594, // I
    113.1594, // L
    128.1741, // K
    131.1926, // M
    147.1766, // F
    97.1167,  // P
    87.0782,  // S
    101.1051, // T
    186.2132, // W
    163.1760, // Y
    99.1326,  // V
];

const WATER_MASS: f64 = 18.01528;

const AA20_HYDRO_KD: [f64; 20] = [
    1.8,  // A
    -4.5, // R
    -3.5, // N
    -3.5, // D
    2.5,  // C
    -3.5, // E
    -3.5, // Q
    -0.4, // G
    -3.2, // H
    4.5,  // I
    3.8,  // L
    -3.9, // K
    1.9,  // M
    2.8,  // F
    -1.6, // P
    -0.8, // S
    -0.7, // T
    -0.9, // W
    -1.3, // Y
    4.2,  // V
];

const PKA_NTERM: f64 = 9.69;
const PKA_CTERM: f64 = 2.34;
const PKA_C: f64 = 8.33;
const PKA_D: f64 = 3.86;
const PKA_E: f64 = 4.25;
const PKA_H: f64 = 6.00;
const PKA_K: f64 = 10.53;
const PKA_R: f64 = 12.48;
const PKA_Y: f64 = 10.07;

#[inline]
fn basic_charge(ph: f64, pka: f64) -> f64 {
    1.0 / (1.0 + 10f64.powf(ph - pka))
}

#[inline]
fn acidic_charge(ph: f64, pka: f64) -> f64 {
    -1.0 / (1.0 + 10f64.powf(pka - ph))
}

#[inline]
fn idx(aa: char) -> usize {
    AA20_INDEX[aa as usize] as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construction_and_reverse() {
        let seq = ProteinSeq::new(b"ACDE".to_vec()).unwrap();
        assert_eq!(seq.as_bytes(), b"ACDE");
        assert_eq!(seq.reverse().as_bytes(), b"EDCA");
        assert!(ProteinSeq::new(b"A#".to_vec()).is_err());
    }

    #[test]
    fn find_and_count() {
        let seq = ProteinSeq::new(b"AAAAA".to_vec()).unwrap();
        assert_eq!(seq.count(b"AA").unwrap(), 2);
        assert_eq!(seq.count_overlap(b"AA").unwrap(), 4);
        assert_eq!(seq.find(b"AAA", 0, 5).unwrap(), Some(0));
        assert_eq!(seq.rfind(b"AAA", 0, 5).unwrap(), Some(2));
        assert!(seq.contains(b"AA").unwrap());
    }

    #[test]
    fn aa_counts_and_entropy() {
        let seq = ProteinSeq::new(b"ACAC".to_vec()).unwrap();
        let counts = seq.aa_counts_20();
        assert_eq!(counts[0], 2); // A
        assert_eq!(counts[4], 2); // C
        let freq = seq.aa_frequencies_20();
        assert!((freq[0] - 0.5).abs() < 1e-12);
        assert!((freq[4] - 0.5).abs() < 1e-12);
        assert!((seq.shannon_entropy() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn molecular_weight() {
        let seq = ProteinSeq::new(b"AC".to_vec()).unwrap();
        let mw = seq.molecular_weight().unwrap();
        assert!((mw - 192.23288).abs() < 1e-6);
    }

    #[test]
    fn hydrophobicity_and_profile() {
        let seq = ProteinSeq::new(b"ACD".to_vec()).unwrap();
        let avg = seq.hydrophobicity().unwrap();
        assert!((avg - 0.2666666667).abs() < 1e-6);
        let profile = seq.hydrophobicity_profile(2).unwrap();
        assert_eq!(profile.len(), 2);
        assert!((profile[0] - 2.15).abs() < 1e-6);
        assert!((profile[1] + 0.5).abs() < 1e-6);
    }

    #[test]
    fn net_charge_and_pi() {
        let seq = ProteinSeq::new(b"AC".to_vec()).unwrap();
        let charge_low = seq.net_charge(2.0).unwrap();
        let charge_high = seq.net_charge(12.0).unwrap();
        assert!(charge_low > 0.0);
        assert!(charge_high < 0.0);
        let pi = seq.isoelectric_point().unwrap();
        let charge_pi = seq.net_charge(pi).unwrap();
        assert!((0.0..=14.0).contains(&pi));
        assert!(charge_pi.abs() < 1e-2);
    }

    #[test]
    fn ambiguity_helpers() {
        let seq = ProteinSeq::new(b"ACBX".to_vec()).unwrap();
        assert!(seq.has_ambiguous());
        assert_eq!(seq.unknown_positions(), vec![2, 3]);
        assert!(seq.validate_strict_20().is_err());
    }
}
