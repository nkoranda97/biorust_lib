from __future__ import annotations

from ._native import (
    AlignmentDNA,
    AlignmentProtein,
    DNARecord,
    DNARecordBatch,
    ProteinRecord,
    ProteinRecordBatch,
)

def msa(
    records: DNARecordBatch
    | list[DNARecord]
    | ProteinRecordBatch
    | list[ProteinRecord],
    *,
    algorithm: str = "clustalo",
    **kwargs,
) -> AlignmentDNA | AlignmentProtein: ...
