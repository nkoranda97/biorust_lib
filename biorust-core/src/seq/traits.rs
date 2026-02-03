use crate::error::BioResult;

pub trait SeqBytes: Clone + Sized {
    fn as_bytes(&self) -> &[u8];
    fn from_bytes(bytes: Vec<u8>) -> BioResult<Self>;
}
