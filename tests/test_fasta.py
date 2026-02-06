import pytest

from biorust import DNA, DNARecord, DNARecordBatch, read_fasta, write_fasta


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


def test_write_fasta_from_records(tmp_path):
    records = [
        DNARecord("seq1", DNA("ACGT"), "some desc"),
        DNARecord("seq2", DNA("A")),
    ]
    path = tmp_path / "out.fasta"
    write_fasta(str(path), records, line_width=2)

    text = path.read_text()
    assert text == ">seq1 some desc\nAC\nGT\n>seq2\nA\n"


def test_write_fasta_from_batch_roundtrip(tmp_path):
    records = [
        DNARecord("seq1", DNA("ACGT"), "some desc"),
        DNARecord("seq2", DNA("A")),
    ]
    batch = DNARecordBatch(records)
    path = tmp_path / "out_batch.fasta"
    write_fasta(str(path), batch)

    roundtrip = read_fasta(str(path))
    assert len(roundtrip) == 2
    assert roundtrip[0].id == "seq1"
    assert roundtrip[0].description == "some desc"
    assert str(roundtrip[0].seq) == "ACGT"
