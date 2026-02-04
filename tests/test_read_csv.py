import pytest

from biorust import (
    DNARecordBatch,
    ProteinRecordBatch,
    csv_columns,
    read_csv,
)


def test_csv_columns(tmp_path):
    data = "id,seq,desc\ns1,ACGT,first\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    assert csv_columns(str(path)) == ["id", "seq", "desc"]


def test_read_csv_dna_by_name(tmp_path):
    data = "id,seq,desc\ns1,ACGT,first\ns2,TT,\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    batch = read_csv(
        str(path), id_col="id", seq_col="seq", desc_col="desc", on_error="raise"
    )
    assert isinstance(batch, DNARecordBatch)
    assert len(batch) == 2
    assert batch[0].id == "s1"
    assert batch[0].description == "first"
    assert str(batch[0].seq) == "ACGT"
    assert batch[1].description is None
    assert batch.skipped == []


def test_read_csv_dna_by_index(tmp_path):
    data = "id,seq\ns1,ACGT\ns2,TT\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    batch = read_csv(str(path), id_col=0, seq_col=1, alphabet="dna", on_error="raise")
    assert isinstance(batch, DNARecordBatch)
    assert [r.id for r in [batch[0], batch[1]]] == ["s1", "s2"]


def test_read_csv_protein(tmp_path):
    data = "id,seq\np1,ACD\np2,GGG\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    batch = read_csv(
        str(path), id_col="id", seq_col="seq", alphabet="protein", on_error="raise"
    )
    assert isinstance(batch, ProteinRecordBatch)
    assert str(batch[0].seq) == "ACD"
    assert batch.skipped == []


def test_read_csv_missing_column(tmp_path):
    data = "id,seq\ns1,ACGT\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    with pytest.raises(ValueError):
        read_csv(str(path), id_col="id", seq_col="missing")


def test_read_csv_invalid_sequence(tmp_path):
    data = "id,seq\ns1,AC#\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    with pytest.raises(ValueError):
        read_csv(str(path), id_col="id", seq_col="seq", on_error="raise")


def test_read_csv_skip_invalid_sequence(tmp_path):
    data = "id,seq\ns1,ACGT\ns2,AC#\ns3,TT\n"
    path = tmp_path / "seqs.csv"
    path.write_text(data)

    batch = read_csv(str(path), id_col="id", seq_col="seq", on_error="skip")
    assert len(batch) == 2
    assert len(batch.skipped) == 1
    skipped = batch.skipped[0]
    assert skipped.row == 2
    assert skipped.id == "s2"
    assert skipped.column
    assert "invalid sequence" in skipped.message
