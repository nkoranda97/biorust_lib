import pytest

from biorust import DNA, DNABatch, Protein, ProteinBatch


def test_dna_batch_index_and_list():
    batch = DNABatch([DNA("AC"), DNA("GT")])

    assert len(batch) == 2
    assert isinstance(batch[0], DNA)
    assert str(batch[1]) == "GT"

    assert [str(s) for s in batch.to_list()] == ["AC", "GT"]
    assert [str(s) for s in list(batch)] == ["AC", "GT"]

    sub = batch[0:2]
    assert isinstance(sub, DNABatch)
    assert [str(s) for s in sub.to_list()] == ["AC", "GT"]


def test_dna_batch_reverse_complements():
    batch = DNABatch([DNA("ATGC"), DNA("AACG")])
    out = batch.reverse_complements()
    assert [str(s) for s in out.to_list()] == ["GCAT", "CGTT"]


def test_dna_batch_reverse_complements_inplace():
    batch = DNABatch([DNA("ATGC"), DNA("AACG")])
    out = batch.reverse_complements(inplace=True)
    assert out is None
    assert [str(s) for s in batch.to_list()] == ["GCAT", "CGTT"]


def test_dna_batch_mutations():
    batch = DNABatch([DNA("A"), DNA("C")])

    batch.append(DNA("G"))
    batch.extend([DNA("T")])
    assert [str(s) for s in batch.to_list()] == ["A", "C", "G", "T"]

    popped = batch.pop()
    assert isinstance(popped, DNA)
    assert str(popped) == "T"

    batch.truncate(2)
    assert [str(s) for s in batch.to_list()] == ["A", "C"]

    batch += DNABatch([DNA("G")])
    assert [str(s) for s in batch.to_list()] == ["A", "C", "G"]

    batch *= 2
    assert [str(s) for s in batch.to_list()] == ["A", "C", "G", "A", "C", "G"]

    batch *= 0
    assert len(batch) == 0

    batch.extend(DNABatch([DNA("AC")]))
    batch.reverse_complements_in_place()
    assert [str(s) for s in batch.to_list()] == ["GT"]

    batch.clear()
    assert len(batch) == 0


def test_protein_batch_index_and_list():
    batch = ProteinBatch([Protein("AC"), Protein("DE")])

    assert len(batch) == 2
    assert isinstance(batch[0], Protein)
    assert str(batch[1]) == "DE"

    assert [str(s) for s in batch.to_list()] == ["AC", "DE"]
    assert [str(s) for s in list(batch)] == ["AC", "DE"]

    sub = batch[0:2]
    assert isinstance(sub, ProteinBatch)
    assert [str(s) for s in sub.to_list()] == ["AC", "DE"]


def test_protein_batch_mutations():
    batch = ProteinBatch([Protein("A"), Protein("C")])

    batch.append(Protein("G"))
    batch.extend([Protein("T")])
    assert [str(s) for s in batch.to_list()] == ["A", "C", "G", "T"]

    popped = batch.pop()
    assert isinstance(popped, Protein)
    assert str(popped) == "T"

    batch.truncate(2)
    assert [str(s) for s in batch.to_list()] == ["A", "C"]

    batch += ProteinBatch([Protein("G")])
    assert [str(s) for s in batch.to_list()] == ["A", "C", "G"]

    batch *= 3
    assert [str(s) for s in batch.to_list()] == [
        "A",
        "C",
        "G",
        "A",
        "C",
        "G",
        "A",
        "C",
        "G",
    ]

    batch *= -1
    assert len(batch) == 0

    batch.clear()
    assert len(batch) == 0


def test_batch_invalid_input():
    with pytest.raises(Exception):
        DNABatch(["AC", DNA("GT")])

    with pytest.raises(Exception):
        ProteinBatch(["AC", Protein("DE")])

    with pytest.raises(Exception):
        batch = DNABatch([DNA("AC")])
        batch.append("AC")

    with pytest.raises(Exception):
        batch = DNABatch([DNA("AC")])
        batch.extend(["AC"])

    with pytest.raises(Exception):
        batch = DNABatch([DNA("AC")])
        batch += ["AC"]

    with pytest.raises(Exception):
        batch = ProteinBatch([Protein("AC")])
        batch.append("AC")

    with pytest.raises(Exception):
        batch = ProteinBatch([Protein("AC")])
        batch.extend(["AC"])

    with pytest.raises(Exception):
        batch = ProteinBatch([Protein("AC")])
        batch += ["AC"]

    with pytest.raises(Exception):
        batch = DNABatch([])
        batch.pop()

    with pytest.raises(Exception):
        batch = ProteinBatch([])
        batch.pop()
