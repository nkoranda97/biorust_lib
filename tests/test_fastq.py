import pytest

from biorust import (
    DNA,
    DNARecord,
    DNARecordBatch,
    ProteinRecordBatch,
    RNARecordBatch,
    read_fastq,
    write_fastq,
)


def test_read_fastq(tmp_path):
    data = "@seq1 some desc\nAC\n+\n!!\n@seq2\nGT\n+\n##\n"
    path = tmp_path / "test.fastq"
    path.write_text(data)

    batch = read_fastq(str(path))
    assert isinstance(batch, DNARecordBatch)
    assert len(batch) == 2
    assert batch[0].id == "seq1"
    assert batch[0].description == "some desc"
    assert str(batch[0].seq) == "AC"
    assert batch[1].id == "seq2"
    assert batch[1].description is None
    assert str(batch[1].seq) == "GT"


def test_read_fastq_invalid_quality_length(tmp_path):
    data = "@seq1\nACGT\n+\n!!!\n"
    path = tmp_path / "bad.fastq"
    path.write_text(data)

    with pytest.raises(ValueError):
        read_fastq(str(path))


def test_write_fastq_from_records(tmp_path):
    records = [
        DNARecord("seq1", DNA("ACGT"), "some desc"),
        DNARecord("seq2", DNA("A")),
    ]
    path = tmp_path / "out.fastq"
    write_fastq(str(path), records, quality_char="J")

    text = path.read_text()
    assert text == "@seq1 some desc\nACGT\n+\nJJJJ\n@seq2\nA\n+\nJ\n"


def test_write_fastq_from_batch_roundtrip(tmp_path):
    records = [
        DNARecord("seq1", DNA("ACGT"), "some desc"),
        DNARecord("seq2", DNA("A")),
    ]
    batch = DNARecordBatch(records)
    path = tmp_path / "out_batch.fastq"
    write_fastq(str(path), batch)

    roundtrip = read_fastq(str(path))
    assert isinstance(roundtrip, DNARecordBatch)
    assert len(roundtrip) == 2
    assert roundtrip[0].id == "seq1"
    assert roundtrip[0].description == "some desc"
    assert str(roundtrip[0].seq) == "ACGT"


def test_read_fastq_auto_detect_rna_and_protein(tmp_path):
    rna_path = tmp_path / "rna.fastq"
    rna_path.write_text("@r1\nACGU\n+\n!!!!\n")
    rna_batch = read_fastq(str(rna_path))
    assert isinstance(rna_batch, RNARecordBatch)

    protein_path = tmp_path / "protein.fastq"
    protein_path.write_text("@p1\nMFW\n+\n###\n")
    protein_batch = read_fastq(str(protein_path))
    assert isinstance(protein_batch, ProteinRecordBatch)
