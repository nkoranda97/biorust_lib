use pyo3::prelude::*;

mod dna;
mod utils;

#[pymodule]
fn _native(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    dna::register(m)?;
    Ok(())
}