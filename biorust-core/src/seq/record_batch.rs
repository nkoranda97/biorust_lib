use crate::error::{BioError, BioResult};
use crate::seq::batch::SeqBatch;
use crate::seq::dna::DnaSeq;
use crate::seq::record::SeqRecord;
use crate::seq::traits::SeqBytes;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecordBatch<S: SeqBytes> {
    ids: Vec<Box<str>>,
    descs: Vec<Option<Box<str>>>,
    seqs: SeqBatch<S>,
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
        Ok(Self {
            ids,
            descs,
            seqs: SeqBatch::new(seqs),
        })
    }

    pub fn from_records(records: Vec<SeqRecord<S>>) -> Self {
        let mut ids = Vec::with_capacity(records.len());
        let mut descs = Vec::with_capacity(records.len());
        let mut seqs = Vec::with_capacity(records.len());

        for record in records {
            ids.push(record.id);
            descs.push(record.desc);
            seqs.push(record.seq);
        }

        Self {
            ids,
            descs,
            seqs: SeqBatch::new(seqs),
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

    pub fn id(&self, i: usize) -> Option<&str> {
        self.ids.get(i).map(|s| s.as_ref())
    }

    pub fn desc(&self, i: usize) -> Option<Option<&str>> {
        self.descs.get(i).map(|d| d.as_deref())
    }

    pub fn seq(&self, i: usize) -> Option<&S> {
        self.seqs.get(i)
    }

    pub fn get_record(&self, i: usize) -> Option<SeqRecordRef<'_, S>> {
        Some(SeqRecordRef {
            id: self.id(i)?,
            desc: self.descs.get(i).and_then(|d| d.as_deref()),
            seq: self.seqs.get(i)?,
        })
    }

    pub fn lengths(&self) -> Vec<usize> {
        self.seqs.lengths()
    }
}

pub struct SeqRecordRef<'a, S: SeqBytes> {
    pub id: &'a str,
    pub desc: Option<&'a str>,
    pub seq: &'a S,
}

impl RecordBatch<DnaSeq> {
    pub fn reverse_complements(&self) -> Self {
        let seqs = self.seqs.reverse_complements().into_vec();
        Self {
            ids: self.ids.clone(),
            descs: self.descs.clone(),
            seqs: SeqBatch::new(seqs),
        }
    }

    pub fn reverse_complements_in_place(&mut self) {
        self.seqs.reverse_complements_in_place();
    }
}
