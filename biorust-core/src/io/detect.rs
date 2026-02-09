/// Sequence-type auto-detection for I/O.
///
/// Rules (deterministic, not probabilistic):
/// - Contains any protein-only character (DEFHIKLMPQRSVWY) → Protein
/// - Contains U but not T → RNA
/// - Otherwise → DNA (the safe default)

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqType {
    Dna,
    Rna,
    Protein,
}

/// Characters that appear in protein sequences but never in DNA/RNA IUPAC alphabets.
const PROTEIN_ONLY: &[u8] = b"DEFHIKLMPQRSVWYdefhiklmpqrsvwy";

/// Detect the sequence type from raw bytes.
///
/// Scans through `bytes` tracking whether any protein-only character is seen,
/// and whether T or U appear. Short-circuits on the first protein-only character.
pub fn detect_seq_type(bytes: &[u8]) -> SeqType {
    let mut has_t = false;
    let mut has_u = false;

    for &b in bytes {
        if PROTEIN_ONLY.contains(&b) {
            return SeqType::Protein;
        }
        match b {
            b'T' | b't' => has_t = true,
            b'U' | b'u' => has_u = true,
            _ => {}
        }
    }

    if has_u && !has_t {
        SeqType::Rna
    } else {
        SeqType::Dna
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_protein() {
        assert_eq!(detect_seq_type(b"MFVFLVLLPLVSS"), SeqType::Protein);
        assert_eq!(detect_seq_type(b"acdefghiklm"), SeqType::Protein);
    }

    #[test]
    fn detect_rna_via_u() {
        assert_eq!(detect_seq_type(b"ACGU"), SeqType::Rna);
        assert_eq!(detect_seq_type(b"acgu"), SeqType::Rna);
        assert_eq!(detect_seq_type(b"AACCGGUU"), SeqType::Rna);
    }

    #[test]
    fn detect_dna_via_t() {
        assert_eq!(detect_seq_type(b"ACGT"), SeqType::Dna);
        assert_eq!(detect_seq_type(b"acgt"), SeqType::Dna);
    }

    #[test]
    fn detect_dna_ambiguous() {
        // Both T and U → defaults to DNA
        assert_eq!(detect_seq_type(b"ACGTU"), SeqType::Dna);
    }

    #[test]
    fn detect_dna_no_tu() {
        // Only A, C, G → defaults to DNA
        assert_eq!(detect_seq_type(b"AACCGG"), SeqType::Dna);
    }

    #[test]
    fn detect_empty() {
        assert_eq!(detect_seq_type(b""), SeqType::Dna);
    }

    #[test]
    fn detect_protein_short_circuits() {
        // First char is protein-only, should short-circuit
        assert_eq!(detect_seq_type(b"MACGT"), SeqType::Protein);
    }

    #[test]
    fn detect_with_iupac_ambiguity_codes() {
        // N, B, X are in both DNA IUPAC and protein — not protein-only
        // Without protein-only chars, this should be DNA
        assert_eq!(detect_seq_type(b"ACGTNB"), SeqType::Dna);
    }
}
