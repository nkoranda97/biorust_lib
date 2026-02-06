use crate::error::{BioError, BioResult};
use std::collections::HashMap;

pub type Qualifiers = HashMap<Box<str>, Vec<Box<str>>>;
pub type Annotations = HashMap<Box<str>, Vec<Box<str>>>;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct FeatureLocation {
    start: usize,
    end: usize,
    strand: Option<i8>,
}

impl FeatureLocation {
    pub fn new(start: usize, end: usize, strand: Option<i8>) -> BioResult<Self> {
        if start > end {
            return Err(BioError::InvalidLocation { start, end });
        }
        if let Some(strand) = strand {
            if strand != -1 && strand != 1 {
                return Err(BioError::InvalidStrand { strand });
            }
        }
        Ok(Self { start, end, strand })
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    pub fn strand(&self) -> Option<i8> {
        self.strand
    }

    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }

    pub fn reverse_complement(&self, len: usize) -> Self {
        debug_assert!(self.end <= len);
        let start = len.saturating_sub(self.end);
        let end = len.saturating_sub(self.start);
        let strand = self.strand.map(|s| -s);
        Self { start, end, strand }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqFeature {
    location: FeatureLocation,
    feature_type: Box<str>,
    qualifiers: Qualifiers,
}

impl SeqFeature {
    pub fn new(feature_type: impl Into<Box<str>>, location: FeatureLocation) -> BioResult<Self> {
        let feature_type = feature_type.into();
        if feature_type.is_empty() {
            return Err(BioError::InvalidFeatureType);
        }
        Ok(Self {
            location,
            feature_type,
            qualifiers: Qualifiers::new(),
        })
    }

    pub fn with_qualifiers(mut self, qualifiers: Qualifiers) -> Self {
        self.qualifiers = qualifiers;
        self
    }

    pub fn location(&self) -> &FeatureLocation {
        &self.location
    }

    pub fn feature_type(&self) -> &str {
        &self.feature_type
    }

    pub fn qualifiers(&self) -> &Qualifiers {
        &self.qualifiers
    }

    pub fn qualifiers_mut(&mut self) -> &mut Qualifiers {
        &mut self.qualifiers
    }

    pub fn set_location(&mut self, location: FeatureLocation) {
        self.location = location;
    }

    pub fn set_feature_type(&mut self, feature_type: impl Into<Box<str>>) -> BioResult<()> {
        let feature_type = feature_type.into();
        if feature_type.is_empty() {
            return Err(BioError::InvalidFeatureType);
        }
        self.feature_type = feature_type;
        Ok(())
    }

    pub fn reverse_complement(&self, len: usize) -> Self {
        let mut out = self.clone();
        out.location = self.location.reverse_complement(len);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_location_validation() {
        assert!(FeatureLocation::new(5, 2, None).is_err());
        assert!(FeatureLocation::new(0, 1, Some(0)).is_err());
        assert!(FeatureLocation::new(0, 1, Some(2)).is_err());
    }

    #[test]
    fn feature_location_basics() {
        let loc = FeatureLocation::new(2, 5, Some(1)).unwrap();
        assert_eq!(loc.start(), 2);
        assert_eq!(loc.end(), 5);
        assert_eq!(loc.strand(), Some(1));
        assert_eq!(loc.len(), 3);
        assert!(!loc.is_empty());
    }

    #[test]
    fn seq_feature_basics() {
        let loc = FeatureLocation::new(0, 3, Some(1)).unwrap();
        let mut quals = Qualifiers::new();
        quals.insert("gene".into(), vec!["abc".into()]);
        let feat = SeqFeature::new("gene", loc)
            .unwrap()
            .with_qualifiers(quals.clone());
        assert_eq!(feat.feature_type(), "gene");
        assert_eq!(feat.qualifiers(), &quals);
    }
}
