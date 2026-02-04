use pyo3::prelude::*;

mod batch;
mod csv;
mod dna;
mod dna_record;
mod dna_record_batch;
mod fasta;
mod protein;
mod protein_record;
mod protein_record_batch;
mod report;
mod seq_shared;
mod utils;

#[pymodule]
fn _native(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    dna::register(m)?;
    dna_record::register(m)?;
    dna_record_batch::register(m)?;
    protein_record::register(m)?;
    protein_record_batch::register(m)?;
    report::register(m)?;
    fasta::register(m)?;
    csv::register(m)?;
    protein::register(m)?;
    batch::register(m)?;
    Ok(())
}
