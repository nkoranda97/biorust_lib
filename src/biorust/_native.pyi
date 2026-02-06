from __future__ import annotations
from typing import Iterable, Iterator, Literal, overload

class Scoring:
    """Scoring parameters for sequence alignment.

    Gap penalties must be <= 0 (negative values). The library uses affine gap penalties:
    score = gap_open + gap_extend * (gap_length - 1)

    Args:
        match_score: Score for matching bases (default: 2)
        mismatch_score: Score for mismatched bases (default: -1)
        gap_open: Penalty for opening a gap (must be <= 0, default: -2.0)
        gap_extend: Penalty for extending a gap (must be <= 0, default: -1.0)
        matrix: Substitution matrix (list of scores or matrix name like 'BLOSUM62', 'EDNAFULL')
        alphabet_size: Required when matrix is a list; size of the alphabet
        end_gap: Whether to use special penalties for terminal gaps
        end_gap_open: Penalty for opening terminal gap (if end_gap=True)
        end_gap_extend: Penalty for extending terminal gap (if end_gap=True)
        use_matrix: Auto-select matrix (EDNAFULL for DNA, BLOSUM62 for protein)
    """
    def __init__(
        self,
        match_score: int | float = ...,
        mismatch_score: int | float = ...,
        gap_open: int | float = ...,
        gap_extend: int | float = ...,
        matrix: list[int] | str | None = ...,
        alphabet_size: int | None = ...,
        end_gap: bool = ...,
        end_gap_open: float | None = ...,
        end_gap_extend: float | None = ...,
        use_matrix: bool = ...,
    ) -> None: ...
    @staticmethod
    def with_matrix(
        matrix: list[int],
        alphabet_size: int,
        gap_open: int | float = ...,
        gap_extend: int | float = ...,
        end_gap: bool = ...,
        end_gap_open: float | None = ...,
        end_gap_extend: float | None = ...,
    ) -> Scoring: ...
    @staticmethod
    def matrix_names() -> list[str]: ...

class AlignmentResult:
    """Result of sequence alignment.

    Contains alignment score, end positions, and optionally CIGAR string and start positions.
    """

    @property
    def score(self) -> float:
        """Alignment score."""
        ...

    @property
    def query_end(self) -> int:
        """End position in query (0-indexed, exclusive)."""
        ...

    @property
    def target_end(self) -> int:
        """End position in target (0-indexed, exclusive)."""
        ...

    @property
    def query_start(self) -> int | None:
        """Start position in query (0-indexed, inclusive). Only available if traceback=True."""
        ...

    @property
    def target_start(self) -> int | None:
        """Start position in target (0-indexed, inclusive). Only available if traceback=True."""
        ...

    @property
    def cigar(self) -> list[tuple[str, int]] | None:
        """CIGAR string as list of (operation, length) tuples. Only available if traceback=True.

        Operations: 'M' (match/mismatch), 'I' (insertion in query), 'D' (deletion from query)
        """
        ...

    def aligned_strings(self) -> tuple[str, str]:
        """Get aligned sequences as strings with gaps represented by '-'.

        Returns:
            Tuple of (aligned_query, aligned_target)

        Raises:
            ValueError: If traceback was not requested
        """
        ...

    def alignment_diagram(self) -> str:
        """Get formatted alignment diagram showing matches and mismatches.

        Returns:
            Multi-line string with query, match indicators, and target

        Raises:
            ValueError: If traceback was not requested
        """
        ...

class DNA:
    """DNA sequence with IUPAC nucleotide alphabet support.

    Supports standard bases (A, C, G, T) and ambiguity codes (N, R, Y, etc.).
    Automatically converts U to T. Case-insensitive on input, preserves case.
    """

    def __init__(self, seq: str | bytes | bytearray | memoryview) -> None:
        """Create a DNA sequence.

        Args:
            seq: Sequence string or bytes-like object

        Raises:
            ValueError: If sequence contains invalid characters
        """
        ...

    def reverse_complement(self) -> DNA:
        """Return the reverse complement of this DNA sequence."""
        ...

    def complement(self) -> DNA:
        """Return the complement of this DNA sequence (without reversing)."""
        ...

    def translate(self) -> Protein:
        """Translate DNA to protein using the standard genetic code.

        Translates in-frame (codons of 3 bases). Partial codons at the end are ignored.
        Ambiguous codons translate to 'X'.
        """
        ...
    def to_bytes(self) -> bytes: ...
    def __len__(self) -> int: ...
    def __bytes__(self) -> bytes: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __hash__(self) -> int: ...
    def __iter__(self) -> Iterator[str]: ...

    # comparisons
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __le__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...

    # concatenation
    def __add__(self, other: DNA) -> DNA: ...

    # mul
    def __mul__(self, other: int) -> "DNA": ...
    def __rmul__(self, other: int) -> "DNA": ...

    # indexing
    @overload
    def __getitem__(self, other: int) -> "DNA": ...
    @overload
    def __getitem__(self, other: slice) -> "DNA": ...
    def __getitem__(self, other): ...
    def count(self, other: str | bytes | bytearray | memoryview | int | DNA) -> int:
        """Count non-overlapping occurrences of a subsequence.

        Args:
            other: Subsequence to count (string, bytes, or DNA object)

        Returns:
            Number of non-overlapping occurrences
        """
        ...

    def count_overlap(
        self, sub: DNA | str | bytes | bytearray | memoryview | int
    ) -> int:
        """Count overlapping occurrences of a subsequence.

        Args:
            sub: Subsequence to count

        Returns:
            Number of overlapping occurrences
        """
        ...

    def __contains__(
        self, other: str | bytes | bytearray | memoryview | int | DNA
    ) -> bool: ...
    def find(
        self,
        sub: str | bytes | bytearray | memoryview | int | "DNA",
        start: int | None = ...,
        end: int | None = ...,
    ) -> int:
        """Find first occurrence of subsequence.

        Args:
            sub: Subsequence to find
            start: Start position (default: 0)
            end: End position (default: sequence length)

        Returns:
            Position of first occurrence, or -1 if not found
        """
        ...
    def index(
        self,
        sub: str | bytes | bytearray | memoryview | int | "DNA",
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def rfind(
        self,
        sub: str | bytes | bytearray | memoryview | int | "DNA",
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def rindex(
        self,
        sub: str | bytes | bytearray | memoryview | int | "DNA",
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def split(
        self,
        sep: str | bytes | bytearray | memoryview | int | "DNA" | None = ...,
        maxsplit: int = ...,
    ) -> list["DNA"]: ...
    def rsplit(
        self,
        sep: str | bytes | bytearray | memoryview | int | "DNA" | None = ...,
        maxsplit: int = ...,
    ) -> list["DNA"]: ...
    def strip(
        self,
        chars: str | bytes | bytearray | memoryview | int | "DNA" | None = ...,
    ) -> "DNA": ...
    def lstrip(
        self,
        chars: str | bytes | bytearray | memoryview | int | "DNA" | None = ...,
    ) -> "DNA": ...
    def rstrip(
        self,
        chars: str | bytes | bytearray | memoryview | int | "DNA" | None = ...,
    ) -> "DNA": ...
    def upper(self) -> "DNA": ...
    def lower(self) -> "DNA": ...
    def startswith(
        self,
        sub: str
        | bytes
        | bytearray
        | memoryview
        | int
        | "DNA"
        | tuple[str | bytes | bytearray | memoryview | int | "DNA", ...],
        start: int | None = ...,
        end: int | None = ...,
    ) -> bool: ...
    def endswith(
        self,
        sub: str
        | bytes
        | bytearray
        | memoryview
        | int
        | "DNA"
        | tuple[str | bytes | bytearray | memoryview | int | "DNA", ...],
        start: int | None = ...,
        end: int | None = ...,
    ) -> bool: ...

class DNARecord:
    def __init__(self, id: str, seq: DNA, desc: str | None = ...) -> None: ...
    @property
    def id(self) -> str: ...
    @property
    def description(self) -> str | None: ...
    @property
    def seq(self) -> DNA: ...
    def __repr__(self) -> str: ...

class DNARecordBatch:
    """Batch of DNA sequence records with IDs and descriptions.

    Efficiently stores multiple DNA sequences along with their metadata.
    Records can be skipped during I/O operations if they contain invalid sequences.
    """

    def __init__(self, records: Iterable[DNARecord] | "DNARecordBatch") -> None:
        """Create a batch from an iterable of DNARecord objects.

        Args:
            records: Iterable of DNARecord objects or another DNARecordBatch
        """
        ...
    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, other: int) -> DNARecord: ...
    @overload
    def __getitem__(self, other: slice) -> "DNARecordBatch": ...
    def __getitem__(self, other): ...
    def ids(self) -> list[str]: ...
    def descriptions(self) -> list[str | None]: ...
    def seqs(self) -> "DNABatch": ...
    @property
    def skipped(self) -> list["SkippedRecord"]: ...
    @overload
    def reverse_complements(
        self, inplace: Literal[False] = ...
    ) -> "DNARecordBatch": ...
    @overload
    def reverse_complements(self, inplace: Literal[True]) -> None: ...

class ProteinRecord:
    def __init__(self, id: str, seq: Protein, desc: str | None = ...) -> None: ...
    @property
    def id(self) -> str: ...
    @property
    def description(self) -> str | None: ...
    @property
    def seq(self) -> Protein: ...
    def __repr__(self) -> str: ...

class ProteinRecordBatch:
    def __init__(
        self, records: Iterable[ProteinRecord] | "ProteinRecordBatch"
    ) -> None: ...
    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, other: int) -> ProteinRecord: ...
    @overload
    def __getitem__(self, other: slice) -> "ProteinRecordBatch": ...
    def __getitem__(self, other): ...
    def ids(self) -> list[str]: ...
    def descriptions(self) -> list[str | None]: ...
    def seqs(self) -> "ProteinBatch": ...
    @property
    def skipped(self) -> list["SkippedRecord"]: ...

class SkippedRecord:
    @property
    def row(self) -> int: ...
    @property
    def id(self) -> str | None: ...
    @property
    def column(self) -> str: ...
    @property
    def message(self) -> str: ...
    def __repr__(self) -> str: ...

class DNABatch:
    """Batch of DNA sequences for efficient processing.

    Provides batch operations and efficient storage of multiple DNA sequences.
    """

    def __init__(self, seqs: Iterable[DNA]) -> None:
        """Create a batch from an iterable of DNA sequences.

        Args:
            seqs: Iterable of DNA objects
        """
        ...
    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, other: int) -> DNA: ...
    @overload
    def __getitem__(self, other: slice) -> "DNABatch": ...
    def __getitem__(self, other): ...
    def __iter__(self): ...
    def to_list(self) -> list[DNA]: ...
    def lengths(self) -> list[int]: ...
    def append(self, seq: DNA) -> None: ...
    def extend(self, seqs: Iterable[DNA] | "DNABatch") -> None: ...
    def clear(self) -> None: ...
    def reserve(self, additional: int) -> None: ...
    def pop(self) -> DNA: ...
    def truncate(self, len: int) -> None: ...
    def __iadd__(self, other: Iterable[DNA] | "DNABatch") -> "DNABatch": ...
    def __imul__(self, n: int) -> "DNABatch": ...
    @overload
    def reverse_complements(self, inplace: Literal[False] = ...) -> "DNABatch": ...
    @overload
    def reverse_complements(self, inplace: Literal[True]) -> None: ...

class Protein:
    """Protein sequence with IUPAC amino acid alphabet support.

    Supports 20 standard amino acids plus B, Z, X, and * (stop).
    Case-insensitive on input.
    """

    def __init__(self, seq: str | bytes | bytearray | memoryview) -> None:
        """Create a protein sequence.

        Args:
            seq: Sequence string or bytes-like object

        Raises:
            ValueError: If sequence contains invalid amino acid codes
        """
        ...

    def to_bytes(self) -> bytes: ...
    def __len__(self) -> int: ...
    def __bytes__(self) -> bytes: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __hash__(self) -> int: ...
    def __iter__(self) -> Iterator[str]: ...
    def reverse(self) -> "Protein":
        """Return the reverse of this protein sequence."""
        ...

    def counts(self) -> list[tuple[str, int]]:
        """Count occurrences of each amino acid.

        Returns:
            List of (amino acid, count) tuples
        """
        ...

    def frequencies(self) -> list[tuple[str, float]]:
        """Calculate frequency of each amino acid.

        Returns:
            List of (amino acid, frequency) tuples
        """
        ...

    def aa_counts_20(self) -> list[tuple[str, int]]:
        """Count occurrences of the 20 standard amino acids only.

        Ambiguous codes (B, Z, X, *) are excluded.

        Returns:
            List of (amino acid, count) tuples in standard order
        """
        ...

    def aa_frequencies_20(self) -> list[tuple[str, float]]:
        """Calculate frequencies of the 20 standard amino acids only.

        Returns:
            List of (amino acid, frequency) tuples
        """
        ...

    def shannon_entropy(self) -> float:
        """Calculate Shannon entropy of the sequence.

        Higher values indicate more diversity in amino acid composition.
        """
        ...

    def molecular_weight(self) -> float:
        """Calculate molecular weight in Daltons.

        Uses average isotopic masses. Includes mass of water added during peptide bond formation.

        Raises:
            ValueError: If sequence contains ambiguous amino acids
        """
        ...

    def hydrophobicity(self) -> float:
        """Calculate mean hydrophobicity using the Kyte-Doolittle scale.

        Returns:
            Mean hydrophobicity value (positive = hydrophobic, negative = hydrophilic)

        Raises:
            ValueError: If sequence contains ambiguous amino acids
        """
        ...

    def hydrophobicity_profile(self, window: int) -> list[float]:
        """Calculate sliding window hydrophobicity profile.

        Args:
            window: Window size for averaging

        Returns:
            List of hydrophobicity values for each window position

        Raises:
            ValueError: If window size is 0 or sequence contains ambiguous amino acids
        """
        ...

    def net_charge(self, ph: float) -> float:
        """Calculate net charge at given pH.

        Uses Henderson-Hasselbalch equation with standard pKa values.

        Args:
            ph: pH value

        Returns:
            Net charge (positive = basic, negative = acidic)

        Raises:
            ValueError: If sequence contains ambiguous amino acids
        """
        ...

    def isoelectric_point(self) -> float:
        """Calculate isoelectric point (pI).

        The pH at which the protein has zero net charge.

        Returns:
            Isoelectric point (pH value)

        Raises:
            ValueError: If sequence contains ambiguous amino acids
        """
        ...
    def validate_strict_20(self) -> None: ...
    def has_ambiguous(self) -> bool: ...
    def unknown_positions(self) -> list[int]: ...

    # comparisons
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __le__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...

    # concatenation
    def __add__(self, other: Protein) -> Protein: ...

    # mul
    def __mul__(self, other: int) -> Protein: ...
    def __rmul__(self, other: int) -> Protein: ...

    # indexing
    @overload
    def __getitem__(self, other: int) -> Protein: ...
    @overload
    def __getitem__(self, other: slice) -> Protein: ...
    def __getitem__(self, other): ...
    def count(
        self, other: str | bytes | bytearray | memoryview | int | Protein
    ) -> int: ...
    def count_overlap(
        self, sub: Protein | str | bytes | bytearray | memoryview | int
    ) -> int: ...
    def __contains__(
        self, other: str | bytes | bytearray | memoryview | int | Protein
    ) -> bool: ...
    def find(
        self,
        sub: str | bytes | bytearray | memoryview | int | Protein,
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def index(
        self,
        sub: str | bytes | bytearray | memoryview | int | Protein,
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def rfind(
        self,
        sub: str | bytes | bytearray | memoryview | int | Protein,
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def rindex(
        self,
        sub: str | bytes | bytearray | memoryview | int | Protein,
        start: int | None = ...,
        end: int | None = ...,
    ) -> int: ...
    def split(
        self,
        sep: str | bytes | bytearray | memoryview | int | Protein | None = ...,
        maxsplit: int = ...,
    ) -> list[Protein]: ...
    def rsplit(
        self,
        sep: str | bytes | bytearray | memoryview | int | Protein | None = ...,
        maxsplit: int = ...,
    ) -> list[Protein]: ...
    def strip(
        self,
        chars: str | bytes | bytearray | memoryview | int | Protein | None = ...,
    ) -> Protein: ...
    def lstrip(
        self,
        chars: str | bytes | bytearray | memoryview | int | Protein | None = ...,
    ) -> Protein: ...
    def rstrip(
        self,
        chars: str | bytes | bytearray | memoryview | int | Protein | None = ...,
    ) -> Protein: ...
    def upper(self) -> Protein: ...
    def lower(self) -> Protein: ...
    def startswith(
        self,
        sub: str
        | bytes
        | bytearray
        | memoryview
        | int
        | Protein
        | tuple[str | bytes | bytearray | memoryview | int | Protein, ...],
        start: int | None = ...,
        end: int | None = ...,
    ) -> bool: ...
    def endswith(
        self,
        sub: str
        | bytes
        | bytearray
        | memoryview
        | int
        | Protein
        | tuple[str | bytes | bytearray | memoryview | int | Protein, ...],
        start: int | None = ...,
        end: int | None = ...,
    ) -> bool: ...

class ProteinBatch:
    def __init__(self, seqs: Iterable[Protein]) -> None: ...
    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, other: int) -> Protein: ...
    @overload
    def __getitem__(self, other: slice) -> "ProteinBatch": ...
    def __getitem__(self, other): ...
    def __iter__(self): ...
    def to_list(self) -> list[Protein]: ...
    def lengths(self) -> list[int]: ...
    def append(self, seq: Protein) -> None: ...
    def extend(self, seqs: Iterable[Protein] | "ProteinBatch") -> None: ...
    def clear(self) -> None: ...
    def reserve(self, additional: int) -> None: ...
    def pop(self) -> Protein: ...
    def truncate(self, len: int) -> None: ...
    def __iadd__(self, other: Iterable[Protein] | "ProteinBatch") -> "ProteinBatch": ...
    def __imul__(self, n: int) -> "ProteinBatch": ...

def complement(seq: str | bytes | bytearray | memoryview | DNA) -> DNA:
    """Compute the complement of a DNA sequence.

    Args:
        seq: DNA sequence (string, bytes, or DNA object)

    Returns:
        Complemented DNA sequence
    """
    ...

def read_fasta(path: str) -> DNARecordBatch:
    """Read DNA sequences from a FASTA file.

    Args:
        path: Path to FASTA file

    Returns:
        DNARecordBatch containing all sequences from the file

    Raises:
        IOError: If file cannot be read
        ValueError: If FASTA format is invalid or sequences contain invalid characters
    """
    ...

@overload
def read_csv(
    path: str,
    *,
    id_col: str | int,
    seq_col: str | int,
    desc_col: str | int | None = ...,
    alphabet: Literal["dna"] = ...,
    on_error: Literal["raise", "skip"] = ...,
) -> DNARecordBatch: ...
@overload
def read_csv(
    path: str,
    *,
    id_col: str | int,
    seq_col: str | int,
    desc_col: str | int | None = ...,
    alphabet: Literal["protein"],
    on_error: Literal["raise", "skip"] = ...,
) -> ProteinRecordBatch: ...
def read_csv(
    path: str,
    *,
    id_col: str | int,
    seq_col: str | int,
    desc_col: str | int | None = ...,
    alphabet: str = ...,
    on_error: str = ...,
) -> DNARecordBatch | ProteinRecordBatch:
    """Read sequences from a CSV file.

    Args:
        path: Path to CSV file
        id_col: Column name or index for sequence IDs
        seq_col: Column name or index for sequences
        desc_col: Optional column name or index for descriptions
        alphabet: Sequence alphabet ('dna' or 'protein', default: 'dna')
        on_error: Error handling ('raise' to fail on invalid sequences, 'skip' to collect in .skipped)

    Returns:
        DNARecordBatch or ProteinRecordBatch containing sequences

    Raises:
        IOError: If file cannot be read
        ValueError: If column not found or on_error='raise' and invalid sequences encountered
    """
    ...

def csv_columns(path: str) -> list[str]: ...
def align_local(
    query: DNA | Protein,
    target: DNA | Protein,
    scoring: Scoring,
    traceback: bool = ...,
) -> AlignmentResult:
    """Perform local (Smith-Waterman) sequence alignment.

    Finds the best local alignment between query and target sequences.
    Uses SIMD acceleration when possible (integer scores, no end gaps).

    Args:
        query: Query sequence (DNA or Protein)
        target: Target sequence (DNA or Protein, must match query type)
        scoring: Scoring parameters
        traceback: If True, compute CIGAR string and start positions (default: False)

    Returns:
        AlignmentResult with score and optionally CIGAR/positions

    Raises:
        ValueError: If sequences are different types or contain invalid characters
    """
    ...

def align_global(
    query: DNA | Protein,
    target: DNA | Protein,
    scoring: Scoring,
    traceback: bool = ...,
) -> AlignmentResult:
    """Perform global (Needleman-Wunsch) sequence alignment.

    Aligns entire sequences end-to-end. Uses SIMD acceleration when possible.

    Args:
        query: Query sequence (DNA or Protein)
        target: Target sequence (DNA or Protein, must match query type)
        scoring: Scoring parameters
        traceback: If True, compute CIGAR string (default: False)

    Returns:
        AlignmentResult with score and optionally CIGAR

    Raises:
        ValueError: If sequences are different types or contain invalid characters
    """
    ...
