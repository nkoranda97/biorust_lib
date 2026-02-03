use crate::error::BioResult;
use crate::seq::dna::ReverseComplement;
use crate::seq::traits::SeqBytes;
use std::ops::Index;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqBatch<S: SeqBytes> {
    seqs: Vec<S>,
}

impl<S: SeqBytes> SeqBatch<S> {
    pub fn new(seqs: Vec<S>) -> Self {
        Self { seqs }
    }

    pub fn from_vec(seqs: Vec<S>) -> Self {
        Self { seqs }
    }

    pub fn as_slice(&self) -> &[S] {
        &self.seqs
    }

    pub fn into_vec(self) -> Vec<S> {
        self.seqs
    }

    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    pub fn get(&self, i: usize) -> Option<&S> {
        self.seqs.get(i)
    }

    pub fn iter(&self) -> impl Iterator<Item = &S> {
        self.seqs.iter()
    }

    pub fn push(&mut self, s: S) {
        self.seqs.push(s);
    }

    pub fn extend<I>(&mut self, iter: I)
    where
        I: IntoIterator<Item = S>,
    {
        self.seqs.extend(iter);
    }

    pub fn clear(&mut self) {
        self.seqs.clear();
    }

    pub fn reserve(&mut self, additional: usize) {
        self.seqs.reserve(additional);
    }

    pub fn pop(&mut self) -> Option<S> {
        self.seqs.pop()
    }

    pub fn truncate(&mut self, len: usize) {
        self.seqs.truncate(len);
    }

    pub fn lengths(&self) -> Vec<usize> {
        self.seqs.iter().map(|seq| seq.as_bytes().len()).collect()
    }

    pub fn to_bytes_vec(&self) -> Vec<Vec<u8>> {
        self.seqs
            .iter()
            .map(|seq| seq.as_bytes().to_vec())
            .collect()
    }

    pub fn map_bytes<F>(&self, f: F) -> BioResult<Self>
    where
        F: Fn(&[u8]) -> Vec<u8>,
    {
        let mut out = Vec::with_capacity(self.seqs.len());
        for seq in &self.seqs {
            let bytes = f(seq.as_bytes());
            out.push(S::from_bytes(bytes)?);
        }
        Ok(Self { seqs: out })
    }

    pub fn map_bytes_in_place<F>(&mut self, f: F) -> BioResult<()>
    where
        F: Fn(&[u8]) -> Vec<u8>,
    {
        let mut out = Vec::with_capacity(self.seqs.len());
        for seq in &self.seqs {
            let bytes = f(seq.as_bytes());
            out.push(S::from_bytes(bytes)?);
        }
        self.seqs = out;
        Ok(())
    }
}

impl<S: SeqBytes> Index<usize> for SeqBatch<S> {
    type Output = S;

    fn index(&self, index: usize) -> &Self::Output {
        &self.seqs[index]
    }
}

impl<S> SeqBatch<S>
where
    S: SeqBytes + ReverseComplement,
{
    pub fn reverse_complements(&self) -> Self {
        let out = self
            .seqs
            .iter()
            .map(|seq| seq.reverse_complement())
            .collect();
        Self { seqs: out }
    }

    pub fn reverse_complements_in_place(&mut self) {
        for seq in &mut self.seqs {
            *seq = seq.reverse_complement();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::DnaSeq;

    #[test]
    fn batch_construction_and_access() {
        let seqs = vec![
            DnaSeq::new(b"AC".to_vec()).unwrap(),
            DnaSeq::new(b"GT".to_vec()).unwrap(),
        ];
        let batch = SeqBatch::new(seqs.clone());

        assert_eq!(batch.len(), 2);
        assert!(!batch.is_empty());
        assert_eq!(batch.as_slice(), seqs.as_slice());
        assert_eq!(batch.get(0), Some(&seqs[0]));

        let collected: Vec<Vec<u8>> = batch.iter().map(|seq| seq.as_bytes().to_vec()).collect();
        assert_eq!(collected, vec![b"AC".to_vec(), b"GT".to_vec()]);
        assert_eq!(batch[1].as_bytes(), b"GT");
    }

    #[test]
    fn batch_lengths_and_bytes() {
        let seqs = vec![
            DnaSeq::new(b"ACG".to_vec()).unwrap(),
            DnaSeq::new(b"T".to_vec()).unwrap(),
        ];
        let batch = SeqBatch::new(seqs);

        assert_eq!(batch.lengths(), vec![3, 1]);
        assert_eq!(batch.to_bytes_vec(), vec![b"ACG".to_vec(), b"T".to_vec()]);
    }

    #[test]
    fn batch_mutations() {
        let mut batch = SeqBatch::new(vec![DnaSeq::new(b"AC".to_vec()).unwrap()]);

        batch.push(DnaSeq::new(b"GT".to_vec()).unwrap());
        assert_eq!(batch.lengths(), vec![2, 2]);

        batch.extend(vec![DnaSeq::new(b"A".to_vec()).unwrap()]);
        assert_eq!(batch.lengths(), vec![2, 2, 1]);

        let popped = batch.pop().unwrap();
        assert_eq!(popped.as_bytes(), b"A");
        assert_eq!(batch.lengths(), vec![2, 2]);

        batch.truncate(1);
        assert_eq!(batch.lengths(), vec![2]);

        let len = batch.len();
        batch.reserve(5);
        assert!(batch.seqs.capacity() >= len + 5);

        batch.clear();
        assert!(batch.is_empty());
    }

    #[test]
    fn batch_map_bytes_success() {
        let seqs = vec![
            DnaSeq::new(b"AC".to_vec()).unwrap(),
            DnaSeq::new(b"GT".to_vec()).unwrap(),
        ];
        let batch = SeqBatch::new(seqs);

        let out = batch
            .map_bytes(|bytes| {
                let mut out = Vec::with_capacity(bytes.len() + 1);
                out.extend_from_slice(bytes);
                out.push(b'A');
                out
            })
            .unwrap();

        assert_eq!(out.to_bytes_vec(), vec![b"ACA".to_vec(), b"GTA".to_vec()]);
    }

    #[test]
    fn batch_map_bytes_in_place_success() {
        let seqs = vec![
            DnaSeq::new(b"AC".to_vec()).unwrap(),
            DnaSeq::new(b"GT".to_vec()).unwrap(),
        ];
        let mut batch = SeqBatch::new(seqs);

        batch
            .map_bytes_in_place(|bytes| {
                let mut out = Vec::with_capacity(bytes.len() + 1);
                out.extend_from_slice(bytes);
                out.push(b'A');
                out
            })
            .unwrap();

        assert_eq!(batch.to_bytes_vec(), vec![b"ACA".to_vec(), b"GTA".to_vec()]);
    }

    #[test]
    fn batch_map_bytes_validation_error() {
        let seqs = vec![DnaSeq::new(b"AC".to_vec()).unwrap()];
        let batch = SeqBatch::new(seqs);

        let err = batch.map_bytes(|_| b"AC#".to_vec());
        assert!(err.is_err());
    }

    #[test]
    fn batch_map_bytes_in_place_validation_error() {
        let seqs = vec![DnaSeq::new(b"AC".to_vec()).unwrap()];
        let mut batch = SeqBatch::new(seqs);

        let err = batch.map_bytes_in_place(|_| b"AC#".to_vec());
        assert!(err.is_err());
        assert_eq!(batch.to_bytes_vec(), vec![b"AC".to_vec()]);
    }

    #[test]
    fn batch_reverse_complements() {
        let seqs = vec![
            DnaSeq::new(b"ATGC".to_vec()).unwrap(),
            DnaSeq::new(b"AACG".to_vec()).unwrap(),
        ];
        let batch = SeqBatch::new(seqs);

        let out = batch.reverse_complements();
        assert_eq!(out.to_bytes_vec(), vec![b"GCAT".to_vec(), b"CGTT".to_vec()]);
    }

    #[test]
    fn batch_reverse_complements_in_place() {
        let seqs = vec![
            DnaSeq::new(b"ATGC".to_vec()).unwrap(),
            DnaSeq::new(b"AACG".to_vec()).unwrap(),
        ];
        let mut batch = SeqBatch::new(seqs);

        batch.reverse_complements_in_place();
        assert_eq!(
            batch.to_bytes_vec(),
            vec![b"GCAT".to_vec(), b"CGTT".to_vec()]
        );
    }
}
