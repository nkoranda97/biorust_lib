import os

import pytest


skip_unless_external = pytest.mark.skipif(
    os.environ.get("BIORUST_RUN_EXTERNAL") != "1",
    reason="set BIORUST_RUN_EXTERNAL=1 to run external tool tests",
)


def test_locator_bogus_path(monkeypatch):
    from biorust._clustalo import find_clustalo

    monkeypatch.setenv("BIORUST_CLUSTALO_PATH", "/tmp/nonexistent_clustalo_binary")
    with pytest.raises(FileNotFoundError, match="does not exist"):
        find_clustalo()


def test_locator_env_override(tmp_path, monkeypatch):
    """Valid path in env var is returned."""
    from biorust._clustalo import find_clustalo

    fake = tmp_path / "clustalo"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)

    monkeypatch.setenv("BIORUST_CLUSTALO_PATH", str(fake))
    result = find_clustalo()
    assert result == fake


@skip_unless_external
def test_msa_smoke():
    from biorust import DNA, DNARecord, DNARecordBatch, msa_clustalo

    records = DNARecordBatch(
        [
            DNARecord("seq1", DNA("ACGTACGT")),
            DNARecord("seq2", DNA("ACGTACGT")),
            DNARecord("seq3", DNA("ACGT")),
        ]
    )

    aln = msa_clustalo(records)

    # ids preserved
    assert set(aln.ids()) == {"seq1", "seq2", "seq3"}

    # all gapped seqs have same length
    assert aln.width > 0
    assert len(aln) == 3

    # only IUPAC + gap chars in output
    valid = set("ACGTRYSWKMBDHVNacgtryswkmbdhvn-.")
    for seq in aln.seqs():
        for ch in str(seq):
            assert ch in valid, f"unexpected char: {ch!r}"


@skip_unless_external
def test_msa_input_list():
    from biorust import DNA, DNARecord, msa_clustalo

    records = [
        DNARecord("a", DNA("ACGTACGT")),
        DNARecord("b", DNA("ACGT")),
    ]

    aln = msa_clustalo(records)
    assert len(aln) == 2
    assert aln.width > 0


@skip_unless_external
def test_msa_protein_smoke():
    from biorust import (
        Protein,
        ProteinRecord,
        ProteinRecordBatch,
        msa_clustalo,
        AlignmentProtein,
    )

    records = ProteinRecordBatch(
        [
            ProteinRecord("seq1", Protein("MKTAYIAKQRQISFVKSH")),
            ProteinRecord("seq2", Protein("MKTAYIAKQRQISFVKSHGGG")),
            ProteinRecord("seq3", Protein("MKTAYIAKQR")),
        ]
    )

    aln = msa_clustalo(records)
    assert isinstance(aln, AlignmentProtein)
    assert set(aln.ids()) == {"seq1", "seq2", "seq3"}
    assert aln.width > 0
    assert len(aln) == 3


@skip_unless_external
def test_msa_generic_dispatch():
    from biorust import DNA, DNARecord, AlignmentDNA, msa

    records = [
        DNARecord("a", DNA("ACGTACGT")),
        DNARecord("b", DNA("ACGT")),
    ]

    aln = msa(records)
    assert isinstance(aln, AlignmentDNA)
    assert len(aln) == 2


@skip_unless_external
def test_msa_generic_dispatch_protein():
    from biorust import Protein, ProteinRecord, AlignmentProtein, msa

    records = [
        ProteinRecord("a", Protein("MKTAYIAKQRQISFVKSH")),
        ProteinRecord("b", Protein("MKTAYIAKQR")),
    ]

    aln = msa(records)
    assert isinstance(aln, AlignmentProtein)
    assert len(aln) == 2


def test_msa_unknown_algorithm():
    from biorust._msa import msa

    with pytest.raises(ValueError, match="unknown algorithm"):
        msa([], algorithm="bogus")


def test_msa_bad_input():
    from biorust._clustalo import _normalize_records

    with pytest.raises(TypeError, match="not str"):
        _normalize_records("ACGT")

    with pytest.raises(TypeError, match="not str"):
        _normalize_records(b"ACGT")
