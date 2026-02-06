# Biorust Codebase Audit — Full Fix Plan

## Context

Broad audit of biorust (Rust bioinformatics library + PyO3 Python bindings) focusing on correctness, type safety, API consistency, and readiness for large datasets. The project aims to be a faster BioPython alternative.

**Scope:** All 3 tiers (correctness, type safety/API, Python ergonomics).
**Gap penalty decision:** Return error on positive gap penalties (explicit, not auto-normalize).

---

## Progress (2026-02-06)

**✅ ALL ITEMS COMPLETE!**

### Completed Items

- ✅ **1.1** — Error messages (DnaSeq/ProteinSeq validation now reports actual char/pos)
- ✅ **1.3** — DnaSeq len/is_empty/reverse (added missing methods)
- ✅ **2.4** — Cigar Display impl (CIGAR string formatting)
- ✅ **1.2** — Scoring Result + validation (replaced assert! panics, added InvalidScoring error)
- ✅ **2.1** — Field encapsulation (EncodedSeq, Cigar, Scoring all have pub(crate) fields + accessors)
- ✅ **2.2 + 2.3** — Remove duplicates, fix in-place semantics (removed from_vec/reverse_complements_inplace alias, fixed RecordBatch bug)
- ✅ **2.5** — Python API consistency (unified DNABatch reverse_complement API)
- ✅ **3.1 + 3.2** — Python iter/hash (added `__iter__` and `__hash__` to DNA and Protein)
- ✅ **3.3** — Docstrings (comprehensive docstrings in _native.pyi for DNA, Protein, Scoring, alignment functions, I/O, and batch types)

### Files Modified

**Rust Core:**

- `biorust-core/src/error.rs` — Added `InvalidScoring` error variant
- `biorust-core/src/seq/dna.rs` — Fixed validation errors, added len/is_empty/reverse, added Hash derive
- `biorust-core/src/seq/protein.rs` — Fixed validation errors, simplified to_string(), added Hash derive
- `biorust-core/src/align/types.rs` — Changed Scoring constructors to return Result, added validation, made fields pub(crate) with accessors, added Cigar Display impl
- `biorust-core/src/align/encode.rs` — Made EncodedSeq fields pub(crate) with accessors
- `biorust-core/src/align/tests.rs` — Updated all Scoring::simple() calls to add .unwrap()
- `biorust-core/src/seq/batch.rs` — Removed from_vec() and reverse_complements_inplace() duplicates
- `biorust-core/src/seq/record_batch.rs` — Fixed reverse_complements_in_place() to actually mutate in place, renamed from _inplace to _in_place

**Python Bindings:**

- `biorust-py/src/align.rs` — Updated to use accessor methods, handle Result from Scoring constructors
- `biorust-py/src/batch.rs` — Removed reverse_complements_in_place() method (kept parameter form only)
- `biorust-py/src/dna.rs` — Added `__hash__`, `__iter__`, and DNAIterator class
- `biorust-py/src/protein.rs` — Added `__hash__`, `__iter__`, and ProteinIterator class
- `biorust-py/src/dna_record_batch.rs` — Updated method call to _in_place
- `src/biorust/_native.pyi` — Added comprehensive docstrings, removed reverse_complements_in_place from DNABatch, added Iterator import

**Tests:**

- `tests/test_batch.py` — Updated to use reverse_complements(inplace=True) instead of reverse_complements_in_place()

### Final Status

**All tests passing:**

- 72 Rust core tests (with simd feature)
- 45 Python tests

**Token usage:** 125,697 / 200,000 (62.8% used)

---

## Tier 1: Correctness Bugs

### 1.1 Fix placeholder error in DnaSeq::new() and ProteinSeq::new()

**Problem:** Both return `BioError::InvalidChar { ch: '?', pos: 0 }` — never reports the actual bad character or position.

**Files:**

- `biorust-core/src/seq/dna.rs:19-24`
- `biorust-core/src/seq/protein.rs:13-17`
- `biorust-core/src/seq/protein.rs:37-41` (ProteinSeq::to_string same issue)

**Change in dna.rs:** Replace:

 ```rust
 pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
     if !dna::iupac_alphabet().is_word(bytes.as_slice()) {
         return Err(BioError::InvalidChar { ch: '?', pos: 0 });
     }
     Ok(Self { bytes })
 }
 ```

**With:**/

```rust
pub fn new(bytes: Vec<u8>) -> BioResult<Self> {
    let alphabet = dna::iupac_alphabet();
    for (pos, &b) in bytes.iter().enumerate() {
        if !alphabet.symbols.contains(b as usize) {
            return Err(BioError::InvalidChar { ch: b as char, pos });
        }
    }
    Ok(Self { bytes })
}
```

Same pattern for `ProteinSeq::new()` using `protein::iupac_alphabet()`.

For `ProteinSeq::to_string()`, replace the `InvalidChar { ch: '?', pos: 0 }` with a more appropriate error or just use `String::from_utf8(self.as_bytes().to_vec()).map_err(...)` — though in practice all valid protein bytes are valid UTF-8 so this error path is unreachable for validated sequences.

**Note:** Check whether `alphabet.symbols` is `pub` — it is (`pub symbols: BitSet` in `biorust-core/src/alphabets/mod.rs:12`).

---

### 1.2 Replace `assert!` panics in Scoring constructors with Result

**Problem:** `Scoring::simple()` and `Scoring::with_matrix()` use `assert!(gap_open <= 0.0)` — panics in library code.

**Files:**

- `biorust-core/src/error.rs` — add new variant
- `biorust-core/src/align/types.rs:70-106` — change constructors
- `biorust-core/src/align/tests.rs` — add `.unwrap()` to all `Scoring::simple()` calls
- `biorust-py/src/align.rs` — propagate Result as PyResult

**Step 1 — Add error variant in error.rs:**

 ```rust
 #[error("invalid scoring parameters: {msg}")]
InvalidScoring { msg: String },
```

**Step 2 — Change `Scoring::simple()` signature from `-> Self` to `-> BioResult<Self>`:**

 ```rust
 pub fn simple(match_score: i16, mismatch_score: i16, gap_open: f32, gap_extend: f32) ->
 BioResult<Self> {
     if gap_open > 0.0 {
         return Err(BioError::InvalidScoring {
             msg: format!("gap_open must be <= 0, got {gap_open}"),
         });
     }
     if gap_extend > 0.0 {
         return Err(BioError::InvalidScoring {
             msg: format!("gap_extend must be <= 0, got {gap_extend}"),
         });
     }
     Ok(Self { ... })
}
```

**Step 3 — Change `Scoring::with_matrix()` similarly, plus validate:**

 ```rust
 if alphabet_size == 0 {
     return Err(BioError::InvalidScoring {
         msg: "alphabet_size must be > 0".into(),
     });
 }
 if matrix.len() != alphabet_size *alphabet_size {
     return Err(BioError::InvalidScoring {
         msg: format!(
             "matrix length {} doesn't match alphabet_size² {}",
             matrix.len(),
             alphabet_size* alphabet_size
         ),
     });
}
```

**Step 4 — Update all call sites.** In `tests.rs`, every `Scoring::simple(...)` becomes `Scoring::simple(...).unwrap()`. In `biorust-py/src/align.rs`, the Scoring `__new__` method already returns `PyResult`, so map the BioError to a Python exception.

---

### 1.3 Add len(), is_empty(), reverse() to DnaSeq

**Problem:** ProteinSeq has these; DnaSeq doesn't.

**File:** `biorust-core/src/seq/dna.rs` — add inside `impl DnaSeq`:

```rust
 pub fn len(&self) -> usize {
     self.bytes.len()
 }

 pub fn is_empty(&self) -> bool {
     self.bytes.is_empty()
 }

 pub fn reverse(&self) -> Self {
     let mut out = self.bytes.clone();
     out.reverse();
        Self { bytes: out }
}
```

---

## Tier 2: Type Safety & API Consistency

### 2.1 Encapsulate public fields on core types

**Problem:** `EncodedSeq`, `Cigar`, and `Scoring` have public fields that allow invalid state.

#### EncodedSeq (`biorust-core/src/align/encode.rs:5-8`)

Change `pub codes` and `pub alphabet_size` to `pub(crate)`. Add:

 ```rust
 pub fn codes(&self) -> &[u8] { &self.codes }
pub fn alphabet_size(&self) -> usize { self.alphabet_size }
```

Internal usages in `align/` module already within the crate, so `pub(crate)` works.

Update `biorust-py/src/align.rs` if it accesses `.codes` or `.alphabet_size` directly — use the accessor methods instead.

#### Cigar (`biorust-core/src/align/types.rs:20-21`)

Change `pub ops` to just `ops`. Add:

 ```rust
 pub fn ops(&self) -> &[(CigarOp, usize)] { &self.ops }
pub fn into_ops(self) -> Vec<(CigarOp, usize)> { self.ops }
```

**Update:**

- `biorust-core/src/align/scalar_ref.rs` — uses `cigar.ops` in `finalize_cigar()`, change to internal construction
- `biorust-py/src/align.rs` — `cigar_to_py()` reads `.ops`, change to `.ops()`

#### Scoring (`biorust-core/src/align/types.rs:58-68`)

Change all 9 `pub` fields to `pub(crate)`. The public API is:

- `Scoring::simple()`, `Scoring::with_matrix()` — constructors
- `Scoring::with_end_gaps()` — builder method
- `score()`, `gap_open_i16()`, `gap_extend_i16()`, `simd_compatible()`, `max_abs_score()` — existing methods

If Python bindings access fields directly, add getter methods as needed.

**Internal access sites to update:**

- `biorust-core/src/align/scalar_ref.rs` — reads `scoring.gap_open`, `scoring.gap_extend`, `scoring.end_gap`, etc.
- `biorust-core/src/align/global_simd.rs` and `local_simd.rs` — reads gap penalties
- `biorust-core/src/align/mod.rs` — reads `scoring.end_gap` for dispatch
- `biorust-py/src/align.rs` — reads fields for `__repr__`

All of these are within the workspace, so `pub(crate)` covers `biorust-core` internal access. For `biorust-py`, add public getter methods where needed:

```rust
 pub fn gap_open(&self) -> f32 { self.gap_open }
 pub fn gap_extend(&self) -> f32 { self.gap_extend }
 pub fn end_gap(&self) -> bool { self.end_gap }
 pub fn end_gap_open(&self) -> f32 { self.end_gap_open }
 pub fn end_gap_extend(&self) -> f32 { self.end_gap_extend }
 pub fn match_score(&self) -> i16 { self.match_score }
pub fn mismatch_score(&self) -> i16 { self.mismatch_score }
pub fn matrix(&self) -> Option<&[i16]> { self.matrix.as_deref() }
```

---

### 2.2 Remove duplicate methods

**File:** `biorust-core/src/seq/batch.rs`

- Remove `from_vec()` (lines 16-18) — identical to `new()`
- Remove `reverse_complements_inplace()` (lines 135-137) — alias of `_in_place`

Check if anything calls them:

- Grep for `from_vec` in the codebase
- Grep for `reverse_complements_inplace` (without underscore before "in")

Update any call sites to use `new()` / `reverse_complements_in_place()`.

---

### 2.3 Fix RecordBatch::reverse_complements_inplace()

**File:** `biorust-core/src/seq/record_batch.rs:113-116`

Replace:

```rust
pub fn reverse_complements_inplace(&mut self) {
    let out = self.seqs.reverse_complements().into_vec();
    self.seqs = SeqBatch::new(out);
}
```

With:

```rust
pub fn reverse_complements_in_place(&mut self) {
    self.seqs.reverse_complements_in_place();
}
```

Also rename from `_inplace` to `_in_place` for consistency.

---

### 2.4 Add Cigar Display impl

**File:** `biorust-core/src/align/types.rs`

```rust
impl std::fmt::Display for Cigar {
     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
         for &(op, len) in &self.ops {
             let ch = match op {
                 CigarOp::Match => 'M',
                 CigarOp::Ins => 'I',
                 CigarOp::Del => 'D',
             };
             write!(f, "{len}{ch}")?;
         }
        Ok(())
    }
}
```

Also expose in Python as `AlignmentResult.cigar_string` property or similar.

---

### 2.5 Unify Python batch reverse_complement API

**Problem:** `DNABatch` has 3 ways (method + parameter + alias), `DNARecordBatch` has 1 (parameter form).

**File:** `biorust-py/src/batch.rs`

Keep only the `reverse_complements(inplace=False)` parameter form. Remove `reverse_complements_in_place()` from Python bindings. Update `src/biorust/_native.pyi` to match.

---

## Tier 3: Python Ergonomics

### 3.1 Add `__iter__` to DNA and Protein classes

**Files:** `biorust-py/src/dna.rs`, `biorust-py/src/protein.rs`

Add `__iter__` that yields single-character DNA/Protein objects. Use PyO3's iterator protocol:

```rust
fn __iter__(slf: PyRef<'_, Self>) -> PyResult<DNAIterator> {
    Ok(DNAIterator {
        bytes: slf.inner.as_bytes().to_vec(),
        index: 0,
    })
}
```

With a separate `#[pyclass]` for the iterator:

```rust
#[pyclass]
struct DNAIterator {
    bytes: Vec<u8>,
    index: usize,
}

#[pymethods]
impl DNAIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> { slf }
    fn __next__(&mut self) -> Option<String> {
        if self.index < self.bytes.len() {
            let ch = self.bytes[self.index] as char;
            self.index += 1;
            Some(ch.to_string())
        } else {
            None
        }
    }
}
```

Register the iterator class in `biorust-py/src/lib.rs`.

Update `src/biorust/_native.pyi` to add `def __iter__(self) -> Iterator[str]: ...`

**Design choice:** Yield `str` (single characters) rather than DNA objects. This matches Python's `str.__iter__()` pattern and avoids creating expensive wrapper objects per base.

---

### 3.2 Add `__hash__` to DNA and Protein classes

**Files:**

- `biorust-core/src/seq/dna.rs` — add Hash derive: `#[derive(Clone, Debug, PartialEq, Eq, Hash)]`
- `biorust-core/src/seq/protein.rs` — same
- `biorust-py/src/dna.rs` — add `__hash__`:

```rust
fn __hash__(&self) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    self.inner.as_bytes().hash(&mut hasher);
    hasher.finish()
}
```

- `biorust-py/src/protein.rs` — same pattern
- `src/biorust/_native.pyi` — add `def __hash__(self) -> int: ...`

---

### 3.3 Add docstrings to .pyi stub

**File:** `src/biorust/_native.pyi`

Add docstrings to the most important types and methods. Priority order:

1. `DNA` class and its key methods (reverse_complement, complement, translate, count, find)
2. `Protein` class and its key methods (molecular_weight, hydrophobicity, isoelectric_point)
3. `Scoring` class (explain gap_open/gap_extend must be <= 0, matrix usage)
4. `align_local`, `align_global` functions
5. `AlignmentResult` properties
6. `read_fasta`, `read_csv` functions
7. Batch types

---

## Implementation Order

1. **1.1** — Error messages (no API breakage, standalone)
2. **1.3** — DnaSeq len/is_empty/reverse (no breakage, standalone)
3. **2.4** — Cigar Display (no breakage, standalone)
4. **1.2** — Scoring Result + validation (breaking change, touches many files)
5. **2.1** — Field encapsulation (breaking, batch all together)
6. **2.2 + 2.3** — Remove duplicates, fix in-place semantics
7. **2.5** — Python API consistency
8. **3.1 + 3.2** — Python iter/hash
9. **3.3** — Docstrings

---

## Verification

After each tier:

- `cargo test -p biorust-core --features simd` — all Rust tests pass
- `cargo clippy --all-targets -- -D warnings` — no warnings
- `uv run maturin develop && uv run pytest -q` — Python tests pass
- `just check` — full pipeline green

---

## Future Work (Not in scope)

- **Streaming I/O for Python** — expose `FastaRecords` iterator to Python, add chunked batch reading
- **Reduce Python batch cloning** — `__iter__`, slicing all clone sequences
- **2-bit DNA encoding** — 4x memory reduction for large genomes
- **Batch alignment API** — align one query against N targets efficiently
- **AlignmentResult field encapsulation** — low priority, output-only struct

---

## Biopython Parity Gaps (Main Missing Features)

- **Scope caveat** — do not include modules that are wrappers around external CLI tools (e.g., MUSCLE/Clustal/EMBOSS/BLAST+ command-line). All APIs must call into our Rust backend only.

| Priority | Area | Scope (Rust backend only) |
| --- | --- | --- |
| P0 | SeqIO write support | FASTA/FASTQ/GenBank/EMBL writers, per-record metadata, batch/streaming writers, line-wrapping controls |
| P0 | RNA sequences | `RNA` type, IUPAC RNA alphabet, complement/rev-comp rules, transcription/back-transcription, RNA translation |
| P0 | SeqRecord features/annotations | `SeqFeature`, `FeatureLocation`, qualifiers, per-record annotations, Python exposure |
| P0 | AlignIO + MSA objects | multiple alignment data model, parse/write Clustal/Stockholm/PHYLIP (at minimum) |
| P1 | PairwiseAligner parity | configurable aligner object (modes, scoring, traceback options, end-gap handling) |
| P1 | Substitution matrix loading | parse standard matrix files (PAM/BLOSUM series), registry of built-in matrices |
| P1 | SeqUtils equivalents | GC/GC skew, melting temp, codon usage, protein analysis helpers not yet covered |
| P1 | SearchIO / BLAST parsers | parse BLAST XML/tabular and HMMER outputs into structured results |
| P1 | Entrez (NCBI E-utilities) | HTTP client with rate limiting, basic response parsing, sequence retrieval |
| P1 | Phylo | tree data model, Newick parsing/writing, traversal utilities |
| P1 | Motifs | PWM/PSSM, consensus, motif file I/O, basic scanning |
| P1 | Restriction enzyme analysis | enzyme catalog, cut site search, digest simulation |
| P1 | Structure I/O (PDB/MMCIF) | structure parsing and basic data model for chains/residues/atoms |
