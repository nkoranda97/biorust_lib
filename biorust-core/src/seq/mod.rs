pub mod batch;
pub mod bytes;
pub mod dna;
pub mod feature;
pub mod gapped_dna;
pub mod gapped_protein;
pub mod protein;
pub mod record;
pub mod record_batch;
pub mod rna;
pub mod traits;

pub use feature::{Annotations, FeatureLocation, Qualifiers, SeqFeature};
pub use record::SeqRecord;
pub use record_batch::{RecordBatch, SeqRecordRef};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TranslationFrame {
    One,
    Two,
    Three,
    Auto,
}

/// Pick the frame (0, 1, or 2) whose translation contains the longest ORF.
/// An ORF is defined as the first M to the next * (or end of sequence).
/// Tiebreak: earliest M start position, then lower frame index.
/// If no M in any frame, returns 0 (frame 1).
pub(crate) fn best_frame_index(proteins: [&[u8]; 3]) -> usize {
    let mut best_frame = 0usize;
    let mut best_len = 0usize;
    let mut best_start = usize::MAX;

    for (frame, protein) in proteins.iter().enumerate() {
        // Find longest ORF: first M to next * or end
        let mut i = 0;
        while i < protein.len() {
            if protein[i] == b'M' {
                let start = i;
                let end = protein[start..].iter().position(|&b| b == b'*');
                let orf_len = match end {
                    Some(e) => e,
                    None => protein.len() - start,
                };
                if orf_len > best_len || (orf_len == best_len && start < best_start) {
                    best_len = orf_len;
                    best_start = start;
                    best_frame = frame;
                }
                break; // only consider first M per frame
            }
            i += 1;
        }
    }

    best_frame
}
