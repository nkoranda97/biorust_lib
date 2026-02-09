"""Tests for translate() and filter_empty() on record batches."""

from biorust import (
    DNA,
    DNARecord,
    DNARecordBatch,
    RNA,
    RNARecord,
    RNARecordBatch,
    Protein,
    ProteinRecord,
    ProteinRecordBatch,
)


# ---- RNA translate ----


def test_rna_record_batch_translate():
    r1 = RNARecord("id1", RNA("AUGAAA"))  # MK
    r2 = RNARecord("id2", RNA("AUGUUU"))  # MF
    batch = RNARecordBatch([r1, r2])

    protein = batch.translate()
    assert isinstance(protein, ProteinRecordBatch)
    assert len(protein) == 2
    assert protein.ids() == ["id1", "id2"]
    seqs = [str(s) for s in protein.seqs().to_list()]
    assert seqs == ["MK", "MF"]


# ---- RNA filter_empty ----


def test_rna_record_batch_filter_empty():
    r1 = RNARecord("id1", RNA("AUGC"))
    r2 = RNARecord("id2", RNA(""))
    r3 = RNARecord("id3", RNA("GCUA"))
    batch = RNARecordBatch([r1, r2, r3])

    filtered = batch.filter_empty()
    assert isinstance(filtered, RNARecordBatch)
    assert len(filtered) == 2
    assert filtered.ids() == ["id1", "id3"]
    assert len(batch) == 3  # original unchanged


def test_rna_record_batch_filter_empty_inplace():
    r1 = RNARecord("id1", RNA("AUGC"))
    r2 = RNARecord("id2", RNA(""))
    batch = RNARecordBatch([r1, r2])

    result = batch.filter_empty(inplace=True)
    assert result is None
    assert len(batch) == 1
    assert batch.ids() == ["id1"]


# ---- Protein filter_empty ----


def test_protein_record_batch_filter_empty():
    r1 = ProteinRecord("id1", Protein("ACDE"))
    r2 = ProteinRecord("id2", Protein(""))
    r3 = ProteinRecord("id3", Protein("FGH"))
    batch = ProteinRecordBatch([r1, r2, r3])

    filtered = batch.filter_empty()
    assert isinstance(filtered, ProteinRecordBatch)
    assert len(filtered) == 2
    assert filtered.ids() == ["id1", "id3"]
    assert len(batch) == 3  # original unchanged


def test_protein_record_batch_filter_empty_inplace():
    r1 = ProteinRecord("id1", Protein("ACDE"))
    r2 = ProteinRecord("id2", Protein(""))
    batch = ProteinRecordBatch([r1, r2])

    result = batch.filter_empty(inplace=True)
    assert result is None
    assert len(batch) == 1
    assert batch.ids() == ["id1"]


# ---- filter_empty with no empties (noop) ----


def test_filter_empty_no_empties():
    r1 = DNARecord("id1", DNA("ATGC"))
    r2 = DNARecord("id2", DNA("GCTA"))
    batch = DNARecordBatch([r1, r2])

    filtered = batch.filter_empty()
    assert len(filtered) == 2


# ---- filter_empty all empties ----


def test_filter_empty_all_empty():
    r1 = DNARecord("id1", DNA(""))
    r2 = DNARecord("id2", DNA(""))
    batch = DNARecordBatch([r1, r2])

    filtered = batch.filter_empty()
    assert len(filtered) == 0
