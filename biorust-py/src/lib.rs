use pyo3::prelude::*;

mod dna;
mod protein;
mod utils;

#[pymodule]
fn _native(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    dna::register(m)?;
    protein::register(m)?;
    Ok(())
}
