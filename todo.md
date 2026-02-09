# To Do

  1. Medium: DistanceMatrix label count is unchecked, which can corrupt indexing or panic if labels don’t match the sequence count.
     Impact: A caller can pass labels.len() != seqs.len() in core APIs, yielding a matrix with an inconsistent n and data length. This can produce
     incorrect distances or out-of-bounds access via get.
     Recommendation: Validate labels.len() == seqs.len() in dna_distance_matrix/protein_distance_matrix (return BioError) or change DistanceMatrix::new
     to check at runtime and return a BioResult.
  2. Medium: Distance computation silently truncates unequal-length sequences because it iterates with zip.
     Evidence: biorust-core/src/phylo/distance.rs:62, biorust-core/src/phylo/distance.rs:94.
     Impact: If sequences are not equal length (especially when calling core APIs directly), trailing characters are ignored and distances are computed
     on the shorter length without warning. This is a correctness risk for core users.
     document “min-length truncation” behavior.
     Test: Add a core test with unequal-length sequences that asserts a clear error.
  3. Medium: Newick serialization does not escape labels, producing invalid output for labels containing whitespace, commas, colons, or parentheses.
     Evidence: biorust-core/src/phylo/newick.rs:20, biorust-core/src/phylo/newick.rs:36.
     Impact: Trees with common IDs like "seq 1" or "A:B" produce malformed Newick, which downstream tools may reject.
     Recommendation: Implement Newick label quoting/escaping (single-quoted labels with escaped quotes) or validate/reject invalid labels at
     construction time.
     Test: Add a test with labels containing spaces/colons and assert correct quoting.
  4. Low (API/UX): Alignment slicing that yields an empty set raises ValueError; Python slicing typically returns an empty alignment (Biopython-style).
     Evidence: biorust-py/src/msa.rs:108, biorust-py/src/msa.rs:302.
     Impact: alignment[0:0] and similar slices error instead of returning an empty alignment, which is surprising for Python users and diverges from
     common expectations.
     Recommendation: Allow empty alignments in slice paths or document the restriction prominently.
     Test: Add Python tests for empty slices on AlignmentDNA and AlignmentProtein.
  5. Low: alignment_diagram treats only '-' as a gap but '.' is also a valid gap for gapped sequences.
     seq/gapped_protein.rs:10.
     Impact: Columns of '.' may be marked as conserved (*) when they should be treated as gaps.
     Recommendation: Treat '.' as a gap in conservation checks.
     Test: Add alignment diagram tests with '.' in aligned sequences.
  6. Low: msa_clustalo allows extra_args to override required flags (-o, --outfmt, etc.), which can break parsing or write output elsewhere.
     Impact: A user can pass extra_args that change output format or path, leading to confusing failures or incorrect parsing.
     Test: Add a unit test that passes a conflicting extra_args and asserts a clean error.
     Evidence: biorust-core/src/phylo/distance.rs:24, biorust-core/src/phylo/distance.rs:42.
     Impact: Misconstructed matrices can panic at runtime.
     Recommendation: Add a runtime check (return error) or make DistanceMatrix::new private and expose a safe constructor.

  Test Gaps

- Newick label escaping/quoting behavior is untested.
- Core phylo distance functions lack tests for unequal sequence lengths and label/sequence length mismatch.
- Alignment slicing empty-result behavior is untested.
- Alignment diagram gaps with '.' are untested.
- msa_clustalo extra_args conflicts are untested.
- MSA integration tests are skipped by default; the suite lacks a lightweight unit-level test for output parsing behavior (without running external
    clustalo).
