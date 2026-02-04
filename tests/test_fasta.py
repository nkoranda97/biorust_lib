import pytest

from biorust import DNARecordBatch, read_fasta


def test_read_fasta(tmp_path):
    data = ">seq1 some desc\nAC\nGT\n>seq2\nA\n"
    path = tmp_path / "test.fasta"
    path.write_text(data)

    batch = read_fasta(str(path))
    assert isinstance(batch, DNARecordBatch)
    assert len(batch) == 2
    assert batch[0].id == "seq1"
    assert batch[0].description == "some desc"
    assert str(batch[0].seq) == "ACGT"
    assert batch[1].id == "seq2"
    assert batch[1].description is None
    assert str(batch[1].seq) == "A"


def test_read_fasta_invalid_char(tmp_path):
    data = ">seq1\nAC#\n"
    path = tmp_path / "bad.fasta"
    path.write_text(data)

    with pytest.raises(ValueError):
        read_fasta(str(path))
