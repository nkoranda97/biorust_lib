use crate::error::BioResult;

pub trait SeqBytes: Clone + Sized {
    fn as_bytes(&self) -> &[u8];
    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self>;

    fn len(&self) -> usize {
        self.as_bytes().len()
    }

    fn is_empty(&self) -> bool {
        self.as_bytes().is_empty()
    }
}
