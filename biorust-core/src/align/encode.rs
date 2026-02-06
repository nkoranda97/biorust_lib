use crate::error::{BioError, BioResult};
use std::sync::LazyLock;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EncodedSeq {
    pub(crate) codes: Vec<u8>,
    pub(crate) alphabet_size: usize,
}

impl EncodedSeq {
    pub fn codes(&self) -> &[u8] {
        &self.codes
    }

    pub fn alphabet_size(&self) -> usize {
        self.alphabet_size
    }

    pub fn len(&self) -> usize {
        self.codes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }
}

const DNA_ALPHABET: &[u8] = b"ATGCSWRYKMBVHDN";
const PROTEIN_ALPHABET: &[u8] = b"ARNDCQEGHILKMFPSTWYVBZX*";

static DNA_MAP: LazyLock<[u8; 256]> = LazyLock::new(|| build_map(DNA_ALPHABET, true));
static PROTEIN_MAP: LazyLock<[u8; 256]> = LazyLock::new(|| build_map(PROTEIN_ALPHABET, true));

fn build_map(alphabet: &[u8], map_lower: bool) -> [u8; 256] {
    let mut map = [255u8; 256];
    for (i, &b) in alphabet.iter().enumerate() {
        map[b as usize] = i as u8;
        if map_lower {
            let lower = b.to_ascii_lowercase();
            map[lower as usize] = i as u8;
        }
    }
    // be forgiving: map U/u to T for DNA
    if alphabet == DNA_ALPHABET {
        map[b'U' as usize] = map[b'T' as usize];
        map[b'u' as usize] = map[b't' as usize];
    }
    map
}

pub fn encode_dna(seq: &[u8]) -> BioResult<EncodedSeq> {
    encode_with_map(seq, &DNA_MAP, DNA_ALPHABET.len())
}

pub fn encode_protein(seq: &[u8]) -> BioResult<EncodedSeq> {
    encode_with_map(seq, &PROTEIN_MAP, PROTEIN_ALPHABET.len())
}

fn encode_with_map(seq: &[u8], map: &[u8; 256], alphabet_size: usize) -> BioResult<EncodedSeq> {
    let mut codes = Vec::with_capacity(seq.len());
    for (pos, &b) in seq.iter().enumerate() {
        let v = map[b as usize];
        if v == 255 {
            return Err(BioError::InvalidChar { ch: b as char, pos });
        }
        codes.push(v);
    }
    Ok(EncodedSeq {
        codes,
        alphabet_size,
    })
}
