# Repository Guidelines

## Project Goals

- Use design of bio-rust as a base for writing a biopython like API
- When appropriate, copy API behavior of biopython
- We want to make a fast biopython while learning some rust

## Project Structure & Module Organization

- `biorust-core/`: Rust library crate; core sequence logic lives in `biorust-core/src/`.
- `biorust-py/`: Rust `pyo3` crate that builds the Python extension; sources in `biorust-py/src/`.
- `src/biorust/`: Python package that re-exports the native module (e.g., `biorust._native`).
- `tests/`: Python tests (`tests/test_*.py`).
- `target/`: Cargo build artifacts (generated).

## Build, Test, and Development Commands

- never run python. Run it through uv as `uv run python`
- `just build` and `just check` for overarching testing.
- `cargo test`: Run Rust tests across the workspace.
- `cargo test -p biorust-core`: Run only the core crate tests.
- `cargo fmt`: Format Rust code with rustfmt.
- `maturin develop`: Build and install the Python extension into the active virtualenv.
- `maturin build`: Build distributable wheels locally.
- `pytest` or `python -m pytest`: Run Python tests in `tests/`.


## Coding Style & Naming Conventions

- Rust: 4-space indentation; `rustfmt` is expected. Use `snake_case` for modules/functions, `PascalCase` for types, and keep public APIs in `biorust-core/src/lib.rs` organized by module.
- Python: 4-space indentation and PEP 8 naming (`snake_case` functions, `PascalCase` classes). The public Python surface should live under `src/biorust/` and wrap `biorust._native` consistently.

## Testing Guidelines

- Rust unit tests live alongside code in `#[cfg(test)]` modules; property-based tests use `proptest` where appropriate. Run with `cargo test`.
- Python tests use `pytest` and should be named `test_*.py`. Add tests for both Rust core behavior and Python bindings when changing APIs.

## Commit & Pull Request Guidelines

- Commit messages in this repo are short, lowercase, and imperative (e.g., `count methods`, `tests`). Keep that style unless a different convention is agreed.
- PRs should include a clear summary, list of tests run, and any relevant issue links. If behavior or Python APIs change, include updated tests and call it out explicitly.

## Configuration Notes

- Python packaging is driven by `pyproject.toml` + `maturin`; Rust workspace config lives in the root `Cargo.toml`.
