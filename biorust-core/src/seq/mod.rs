pub mod batch;
pub mod bytes;
pub mod dna;
pub mod protein;
pub mod record;
pub mod record_batch;
pub mod traits;

pub use record::SeqRecord;
pub use record_batch::{RecordBatch, SeqRecordRef};
