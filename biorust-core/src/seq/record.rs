use crate::seq::feature::{Annotations, SeqFeature};
use crate::seq::traits::SeqBytes;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqRecord<S: SeqBytes> {
    pub id: Box<str>,
    pub desc: Option<Box<str>>,
    pub seq: S,
    pub features: Vec<SeqFeature>,
    pub annotations: Annotations,
}

impl<S: SeqBytes> SeqRecord<S> {
    pub fn new(id: impl Into<Box<str>>, seq: S) -> Self {
        Self {
            id: id.into(),
            desc: None,
            seq,
            features: Vec::new(),
            annotations: Annotations::new(),
        }
    }

    pub fn with_desc(mut self, desc: impl Into<Box<str>>) -> Self {
        self.desc = Some(desc.into());
        self
    }

    pub fn with_features(mut self, features: Vec<SeqFeature>) -> Self {
        self.features = features;
        self
    }

    pub fn with_annotations(mut self, annotations: Annotations) -> Self {
        self.annotations = annotations;
        self
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn desc(&self) -> Option<&str> {
        self.desc.as_deref()
    }

    pub fn seq(&self) -> &S {
        &self.seq
    }

    pub fn features(&self) -> &[SeqFeature] {
        &self.features
    }

    pub fn features_mut(&mut self) -> &mut Vec<SeqFeature> {
        &mut self.features
    }

    pub fn annotations(&self) -> &Annotations {
        &self.annotations
    }

    pub fn annotations_mut(&mut self) -> &mut Annotations {
        &mut self.annotations
    }

    pub fn into_seq(self) -> S {
        self.seq
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::dna::DnaSeq;
    use crate::seq::feature::{Annotations, FeatureLocation, SeqFeature};

    #[test]
    fn record_features_annotations_roundtrip() {
        let seq = DnaSeq::new(b"ATGC".to_vec()).unwrap();
        let loc = FeatureLocation::new(0, 2, Some(1)).unwrap();
        let feature = SeqFeature::new("gene", loc).unwrap();
        let mut ann = Annotations::new();
        ann.insert("source".into(), vec!["test".into()]);

        let record = SeqRecord::new("id1", seq)
            .with_features(vec![feature.clone()])
            .with_annotations(ann.clone());

        assert_eq!(record.features(), &[feature]);
        assert_eq!(record.annotations(), &ann);
    }
}
