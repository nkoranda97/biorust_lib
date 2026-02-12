use crate::error::{BioError, BioResult};
use crate::seq::batch::SeqBatch;
use crate::seq::dna::DnaSeq;
use crate::seq::feature::{Annotations, SeqFeature};
use crate::seq::protein::ProteinSeq;
use crate::seq::record::SeqRecord;
use crate::seq::rna::RnaSeq;
use crate::seq::traits::SeqBytes;
use crate::seq::TranslationFrame;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecordBatch<S: SeqBytes> {
    ids: Vec<Box<str>>,
    descs: Vec<Option<Box<str>>>,
    seqs: SeqBatch<S>,
    features: Vec<Vec<SeqFeature>>,
    annotations: Vec<Annotations>,
}

impl<S: SeqBytes> RecordBatch<S> {
    pub fn new(ids: Vec<Box<str>>, descs: Vec<Option<Box<str>>>, seqs: Vec<S>) -> BioResult<Self> {
        if ids.len() != seqs.len() || descs.len() != seqs.len() {
            return Err(BioError::RecordBatchLenMismatch {
                ids: ids.len(),
                descs: descs.len(),
                seqs: seqs.len(),
            });
        }
        let mut features = Vec::with_capacity(seqs.len());
        let mut annotations = Vec::with_capacity(seqs.len());
        for _ in 0..seqs.len() {
            features.push(Vec::new());
            annotations.push(Annotations::new());
        }
        Ok(Self {
            ids,
            descs,
            seqs: SeqBatch::new(seqs),
            features,
            annotations,
        })
    }

    pub fn new_with_meta(
        ids: Vec<Box<str>>,
        descs: Vec<Option<Box<str>>>,
        seqs: Vec<S>,
        features: Vec<Vec<SeqFeature>>,
        annotations: Vec<Annotations>,
    ) -> BioResult<Self> {
        if ids.len() != seqs.len()
            || descs.len() != seqs.len()
            || features.len() != seqs.len()
            || annotations.len() != seqs.len()
        {
            return Err(BioError::RecordBatchLenMismatch {
                ids: ids.len(),
                descs: descs.len(),
                seqs: seqs.len(),
            });
        }
        Ok(Self {
            ids,
            descs,
            seqs: SeqBatch::new(seqs),
            features,
            annotations,
        })
    }

    pub fn from_records(records: Vec<SeqRecord<S>>) -> Self {
        let mut ids = Vec::with_capacity(records.len());
        let mut descs = Vec::with_capacity(records.len());
        let mut seqs = Vec::with_capacity(records.len());
        let mut features = Vec::with_capacity(records.len());
        let mut annotations = Vec::with_capacity(records.len());

        for record in records {
            ids.push(record.id);
            descs.push(record.desc);
            seqs.push(record.seq);
            features.push(record.features);
            annotations.push(record.annotations);
        }

        Self {
            ids,
            descs,
            seqs: SeqBatch::new(seqs),
            features,
            annotations,
        }
    }

    pub fn len(&self) -> usize {
        self.seqs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seqs.is_empty()
    }

    pub fn ids(&self) -> &[Box<str>] {
        &self.ids
    }

    pub fn descs(&self) -> &[Option<Box<str>>] {
        &self.descs
    }

    pub fn seqs(&self) -> &SeqBatch<S> {
        &self.seqs
    }

    pub fn seqs_mut(&mut self) -> &mut SeqBatch<S> {
        &mut self.seqs
    }

    pub fn features(&self) -> &[Vec<SeqFeature>] {
        &self.features
    }

    pub fn annotations(&self) -> &[Annotations] {
        &self.annotations
    }

    pub fn id(&self, i: usize) -> Option<&str> {
        self.ids.get(i).map(|s| s.as_ref())
    }

    pub fn desc(&self, i: usize) -> Option<Option<&str>> {
        self.descs.get(i).map(|d| d.as_deref())
    }

    pub fn seq(&self, i: usize) -> Option<&S> {
        self.seqs.get(i)
    }

    pub fn features_at(&self, i: usize) -> Option<&[SeqFeature]> {
        self.features.get(i).map(|f| f.as_slice())
    }

    pub fn annotations_at(&self, i: usize) -> Option<&Annotations> {
        self.annotations.get(i)
    }

    pub fn get_record(&self, i: usize) -> Option<SeqRecordRef<'_, S>> {
        Some(SeqRecordRef {
            id: self.id(i)?,
            desc: self.descs.get(i).and_then(|d| d.as_deref()),
            seq: self.seqs.get(i)?,
            features: self.features_at(i)?,
            annotations: self.annotations_at(i)?,
        })
    }

    pub fn lengths(&self) -> Vec<usize> {
        self.seqs.lengths()
    }

    /// Return a new batch containing only records whose sequence is non-empty.
    pub fn filter_empty(&self) -> Self {
        let mut ids = Vec::new();
        let mut descs = Vec::new();
        let mut seqs = Vec::new();
        let mut features = Vec::new();
        let mut annotations = Vec::new();

        for (i, seq) in self.seqs.as_slice().iter().enumerate() {
            if !seq.as_bytes().is_empty() {
                ids.push(self.ids[i].clone());
                descs.push(self.descs[i].clone());
                seqs.push(seq.clone());
                features.push(self.features[i].clone());
                annotations.push(self.annotations[i].clone());
            }
        }

        Self {
            ids,
            descs,
            seqs: SeqBatch::new(seqs),
            features,
            annotations,
        }
    }

    /// Remove records with empty sequences in place.
    pub fn filter_empty_in_place(&mut self) {
        let keep: Vec<bool> = self
            .seqs
            .as_slice()
            .iter()
            .map(|s| !s.as_bytes().is_empty())
            .collect();

        fn retain_by_mask<T>(v: &mut Vec<T>, keep: &[bool]) {
            let mut iter = keep.iter();
            v.retain(|_| *iter.next().unwrap());
        }

        retain_by_mask(&mut self.ids, &keep);
        retain_by_mask(&mut self.descs, &keep);
        retain_by_mask(&mut self.features, &keep);
        retain_by_mask(&mut self.annotations, &keep);

        let mut seqs = self.seqs.as_slice().to_vec();
        retain_by_mask(&mut seqs, &keep);
        self.seqs = SeqBatch::new(seqs);
    }
}

pub struct SeqRecordRef<'a, S: SeqBytes> {
    pub id: &'a str,
    pub desc: Option<&'a str>,
    pub seq: &'a S,
    pub features: &'a [SeqFeature],
    pub annotations: &'a Annotations,
}

impl RecordBatch<DnaSeq> {
    /// Translate all DNA sequences to protein, preserving IDs, descriptions,
    /// and annotations. Features are cleared since nucleotide coordinates
    /// don't map to protein coordinates.
    pub fn translate(&self) -> BioResult<RecordBatch<ProteinSeq>> {
        let seqs = self.seqs.translate()?.into_vec();
        let empty_features: Vec<Vec<SeqFeature>> = vec![Vec::new(); seqs.len()];
        Ok(RecordBatch {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
            features: empty_features,
            annotations: self.annotations.clone(),
        })
    }

    pub fn translate_frame(&self, frame: TranslationFrame) -> BioResult<RecordBatch<ProteinSeq>> {
        let seqs = self.seqs.translate_frame(frame)?.into_vec();
        let empty_features: Vec<Vec<SeqFeature>> = vec![Vec::new(); seqs.len()];
        Ok(RecordBatch {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
            features: empty_features,
            annotations: self.annotations.clone(),
        })
    }

    pub fn reverse_complements(&self) -> Self {
        let seqs = self.seqs.reverse_complements().into_vec();
        let mut features = Vec::with_capacity(self.features.len());
        for (idx, feats) in self.features.iter().enumerate() {
            let len = self.seqs.as_slice()[idx].as_bytes().len();
            let out = feats
                .iter()
                .map(|feat| feat.reverse_complement(len))
                .collect();
            features.push(out);
        }
        Self {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
            features,
            annotations: self.annotations.clone(),
        }
    }

    pub fn reverse_complements_in_place(&mut self) {
        for (idx, feats) in self.features.iter_mut().enumerate() {
            let len = self.seqs.as_slice()[idx].as_bytes().len();
            for feat in feats.iter_mut() {
                *feat = feat.reverse_complement(len);
            }
        }
        self.seqs.reverse_complements_in_place();
    }
}

impl RecordBatch<RnaSeq> {
    /// Translate all RNA sequences to protein, preserving IDs, descriptions,
    /// and annotations. Features are cleared since nucleotide coordinates
    /// don't map to protein coordinates.
    pub fn translate(&self) -> BioResult<RecordBatch<ProteinSeq>> {
        let seqs = self.seqs.translate()?.into_vec();
        let empty_features: Vec<Vec<SeqFeature>> = vec![Vec::new(); seqs.len()];
        Ok(RecordBatch {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
            features: empty_features,
            annotations: self.annotations.clone(),
        })
    }

    pub fn translate_frame(&self, frame: TranslationFrame) -> BioResult<RecordBatch<ProteinSeq>> {
        let seqs = self.seqs.translate_frame(frame)?.into_vec();
        let empty_features: Vec<Vec<SeqFeature>> = vec![Vec::new(); seqs.len()];
        Ok(RecordBatch {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
            features: empty_features,
            annotations: self.annotations.clone(),
        })
    }

    pub fn reverse_complements(&self) -> Self {
        let seqs = self.seqs.reverse_complements().into_vec();
        let mut features = Vec::with_capacity(self.features.len());
        for (idx, feats) in self.features.iter().enumerate() {
            let len = self.seqs.as_slice()[idx].as_bytes().len();
            let out = feats
                .iter()
                .map(|feat| feat.reverse_complement(len))
                .collect();
            features.push(out);
        }
        Self {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
            features,
            annotations: self.annotations.clone(),
        }
    }

    pub fn reverse_complements_in_place(&mut self) {
        for (idx, feats) in self.features.iter_mut().enumerate() {
            let len = self.seqs.as_slice()[idx].as_bytes().len();
            for feat in feats.iter_mut() {
                *feat = feat.reverse_complement(len);
            }
        }
        self.seqs.reverse_complements_in_place();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::DnaSeq;
    use crate::seq::feature::{Annotations, FeatureLocation, SeqFeature};
    use crate::seq::record::SeqRecord;

    #[test]
    fn batch_preserves_features_annotations() {
        let seq = DnaSeq::new(b"ATGC".to_vec()).unwrap();
        let loc = FeatureLocation::new(0, 2, Some(1)).unwrap();
        let feature = SeqFeature::new("gene", loc).unwrap();
        let mut ann = Annotations::new();
        ann.insert("source".into(), vec!["test".into()]);

        let record = SeqRecord::new("id1", seq)
            .with_features(vec![feature.clone()])
            .with_annotations(ann.clone());

        let batch = RecordBatch::from_records(vec![record]);
        assert_eq!(batch.features_at(0).unwrap(), &[feature.clone()]);
        assert_eq!(batch.annotations_at(0).unwrap(), &ann);

        let record_ref = batch.get_record(0).unwrap();
        assert_eq!(record_ref.features, &[feature]);
        assert_eq!(record_ref.annotations, &ann);
    }

    #[test]
    fn translate_dna_batch() {
        let seq1 = DnaSeq::new(b"ATGAAA".to_vec()).unwrap(); // MK
        let seq2 = DnaSeq::new(b"ATGTTT".to_vec()).unwrap(); // MF
        let r1 = SeqRecord::new("id1", seq1);
        let r2 = SeqRecord::new("id2", seq2);
        let batch = RecordBatch::from_records(vec![r1, r2]);

        let protein_batch = batch.translate().unwrap();
        assert_eq!(protein_batch.len(), 2);
        assert_eq!(protein_batch.ids(), batch.ids());
        assert_eq!(protein_batch.seq(0).unwrap().as_bytes(), b"MK");
        assert_eq!(protein_batch.seq(1).unwrap().as_bytes(), b"MF");
    }

    #[test]
    fn filter_empty_records() {
        let seq1 = DnaSeq::new(b"ATGC".to_vec()).unwrap();
        let seq2 = DnaSeq::new(Vec::new()).unwrap();
        let seq3 = DnaSeq::new(b"GCTA".to_vec()).unwrap();
        let r1 = SeqRecord::new("id1", seq1);
        let r2 = SeqRecord::new("id2", seq2);
        let r3 = SeqRecord::new("id3", seq3);
        let batch = RecordBatch::from_records(vec![r1, r2, r3]);

        let filtered = batch.filter_empty();
        assert_eq!(filtered.len(), 2);
        assert_eq!(filtered.id(0).unwrap(), "id1");
        assert_eq!(filtered.id(1).unwrap(), "id3");
    }

    #[test]
    fn filter_empty_in_place_records() {
        let seq1 = DnaSeq::new(b"ATGC".to_vec()).unwrap();
        let seq2 = DnaSeq::new(Vec::new()).unwrap();
        let seq3 = DnaSeq::new(b"GCTA".to_vec()).unwrap();
        let r1 = SeqRecord::new("id1", seq1);
        let r2 = SeqRecord::new("id2", seq2);
        let r3 = SeqRecord::new("id3", seq3);
        let mut batch = RecordBatch::from_records(vec![r1, r2, r3]);

        batch.filter_empty_in_place();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.id(0).unwrap(), "id1");
        assert_eq!(batch.id(1).unwrap(), "id3");
    }

    #[test]
    fn reverse_complement_updates_features() {
        let seq = DnaSeq::new(b"ATGC".to_vec()).unwrap();
        let loc = FeatureLocation::new(0, 2, Some(1)).unwrap();
        let feature = SeqFeature::new("gene", loc).unwrap();
        let record = SeqRecord::new("id1", seq).with_features(vec![feature]);
        let batch = RecordBatch::from_records(vec![record]);

        let rc = batch.reverse_complements();
        let rc_feature = &rc.features()[0][0];
        assert_eq!(rc_feature.location().start(), 2);
        assert_eq!(rc_feature.location().end(), 4);
        assert_eq!(rc_feature.location().strand(), Some(-1));
    }
}
