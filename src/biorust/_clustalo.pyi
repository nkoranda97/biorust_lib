from __future__ import annotations

from pathlib import Path

from ._native import (
    AlignmentDNA,
    AlignmentProtein,
    DNARecord,
    DNARecordBatch,
    ProteinRecord,
    ProteinRecordBatch,
)

def find_clustalo() -> Path: ...
def msa_clustalo(
    records: DNARecordBatch
    | list[DNARecord]
    | ProteinRecordBatch
    | list[ProteinRecord],
    *,
    threads: int | None = None,
    extra_args: list[str] | None = None,
) -> AlignmentDNA | AlignmentProtein: ...
