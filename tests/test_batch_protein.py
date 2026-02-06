import pytest

from biorust import Protein, ProteinBatch


def test_protein_batch_copy_slice_take_concat():
    batch = ProteinBatch([Protein("A"), Protein("CC"), Protein("GGG")])

    copied = batch.copy()
    assert copied is not batch
    assert [str(s) for s in copied.to_list()] == ["A", "CC", "GGG"]

    sliced = batch.slice(0, 3, 2)
    assert [str(s) for s in sliced.to_list()] == ["A", "GGG"]

    sliced_neg = batch.slice(-2, None, 1)
    assert [str(s) for s in sliced_neg.to_list()] == ["CC", "GGG"]

    taken = batch.take([1, -1])
    assert [str(s) for s in taken.to_list()] == ["CC", "GGG"]

    concat = batch.concat()
    assert str(concat) == "ACCGGG"

    with pytest.raises(ValueError):
        ProteinBatch([]).concat()


def test_protein_batch_filter_by_len():
    batch = ProteinBatch([Protein("A"), Protein("CC"), Protein("GGG")])

    filtered = batch.filter_by_len(2, None)
    assert [str(s) for s in filtered.to_list()] == ["CC", "GGG"]

    out = batch.filter_by_len(2, None, inplace=True)
    assert out is None
    assert [str(s) for s in batch.to_list()] == ["CC", "GGG"]


def test_protein_batch_count_contains():
    batch = ProteinBatch([Protein("AAC"), Protein("TT")])

    assert batch.count(Protein("A")) == [2, 0]
    assert batch.contains(Protein("TT")) == [False, True]

    with pytest.raises(TypeError):
        batch.count("A")

    with pytest.raises(TypeError):
        batch.contains("TT")


def test_protein_batch_take_and_slice_errors():
    batch = ProteinBatch([Protein("A")])

    with pytest.raises(IndexError):
        batch.take([3])

    with pytest.raises(ValueError):
        batch.slice(step=0)
