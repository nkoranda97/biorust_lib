# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Goal

A fast BioPython alternative built in Rust with Python bindings via PyO3. Copies BioPython API behavior where appropriate.

## Build & Test Commands

```bash
just check          # Full pipeline: fmt, lint, test (Rust + Python), build Python extension
just dev            # Quick: build and install Python extension into virtualenv
cargo test --all    # All Rust tests (simd feature is on by default)
cargo test -p biorust-core              # Core crate only
cargo test -p biorust-core -- test_name # Single Rust test
uv run maturin develop                  # Build Python extension
uv run pytest -q                        # Python tests
uv run pytest tests/test_dna.py -q      # Single Python test file
```

Always run Python through `uv run` (not bare `python` or `pytest`).

Linting: `cargo fmt --all`, `cargo clippy --all-targets -- -D warnings`, `uvx ruff format .`

## Architecture

**Workspace:** Two Rust crates + a Python package.

- `biorust-core/` — Core Rust library. Five modules: `align`, `alphabets`, `seq`, `io`, `error`
- `biorust-py/` — PyO3 bindings wrapping core types into Python classes
- `src/biorust/` — Python package re-exporting `biorust._native`
- `tests/` — Python pytest tests

**Sequence types:** `DnaSeq` and `ProteinSeq` are separate concrete types (not generic over alphabet). Both implement the `SeqBytes` trait. Sequences store raw bytes internally; alignment uses `EncodedSeq` with integer-coded symbols.

**Batches:** `SeqBatch<S>` (Vec wrapper) and `RecordBatch<S>` (parallel vecs of ids, descs, sequences) provide batch operations over `SeqRecord<S>`.

**Alignment (`biorust-core/src/align/`):**

- `mod.rs` dispatches to SIMD or scalar based on `traceback` flag and `simd` feature
- SIMD (score-only): 16-lane `i16x16` via `wide` crate, profile-based DP. `simd_safe_len()` guards against i16 overflow
- Scalar (with traceback): Full affine gap with bit-packed u8 direction array
- Gap penalty convention: SIMD uses positive penalties (subtracted); scalar normalizes to negative internally. Users may pass either sign (BLAST convention supported)
- CIGAR does NOT merge consecutive same-op blocks — preserves gap-open boundaries for pathological cases where |gap_open| < |gap_extend|

**I/O:** FASTA and CSV readers return `RecordBatch<S>`. Error handling via `OnError::Raise` (fail fast) or `OnError::Skip` (collect `SkippedRecord`s in `ReadReport`).

## Feature Flags

- `simd` (default: enabled) — Gates SIMD alignment implementations. Without it, all alignment uses scalar reference code.

## Coding Conventions

- Rust: `rustfmt` with max_width=100, 4-space indent, edition 2021
- Python: PEP 8, `ruff` formatted
- Commit messages: short, lowercase, imperative (e.g., `count methods`, `fasta io`)
- Rust tests in `#[cfg(test)]` modules alongside code; property tests with `proptest`
- Python tests named `test_*.py` in `tests/`
