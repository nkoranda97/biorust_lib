use pyo3::prelude::*;

mod batch;
mod dna;
mod protein;
mod seq_shared;
mod utils;

#[pymodule]
fn _native(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    dna::register(m)?;
    protein::register(m)?;
    batch::register(m)?;
    Ok(())
}
