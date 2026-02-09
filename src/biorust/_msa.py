"""Generic MSA dispatch."""

from __future__ import annotations


def msa(records, *, algorithm="clustalo", **kwargs):
    """Run multiple sequence alignment.

    Parameters
    ----------
    records : DNARecordBatch | list[DNARecord] | ProteinRecordBatch | list[ProteinRecord]
        Input sequences to align.
    algorithm : str
        Alignment algorithm (default: ``"clustalo"``).
    **kwargs
        Forwarded to the algorithm-specific function.

    Returns
    -------
    AlignmentDNA | AlignmentProtein
        Aligned (gapped) sequences.
    """
    if algorithm == "clustalo":
        from biorust._clustalo import msa_clustalo

        return msa_clustalo(records, **kwargs)

    raise ValueError(f"unknown algorithm {algorithm!r}, supported: 'clustalo'")
