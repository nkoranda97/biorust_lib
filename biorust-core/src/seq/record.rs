use crate::seq::traits::SeqBytes;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqRecord<S: SeqBytes> {
    pub id: Box<str>,
    pub desc: Option<Box<str>>,
    pub seq: S,
}

impl<S: SeqBytes> SeqRecord<S> {
    pub fn new(id: impl Into<Box<str>>, seq: S) -> Self {
        Self {
            id: id.into(),
            desc: None,
            seq,
        }
    }

    pub fn with_desc(mut self, desc: impl Into<Box<str>>) -> Self {
        self.desc = Some(desc.into());
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

    pub fn into_seq(self) -> S {
        self.seq
    }
}
