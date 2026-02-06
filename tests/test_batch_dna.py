import pytest

from biorust import DNA, DNABatch, ProteinBatch, RNABatch


def test_dna_batch_copy_slice_take_concat():
    batch = DNABatch([DNA("A"), DNA("CC"), DNA("GGG"), DNA("TTTT")])

    copied = batch.copy()
    assert copied is not batch
    assert [str(s) for s in copied.to_list()] == ["A", "CC", "GGG", "TTTT"]

    sliced = batch.slice(1, None, 2)
    assert [str(s) for s in sliced.to_list()] == ["CC", "TTTT"]

    sliced_neg = batch.slice(-3, -1, 1)
    assert [str(s) for s in sliced_neg.to_list()] == ["CC", "GGG"]

    taken = batch.take([2, 0, -1])
    assert [str(s) for s in taken.to_list()] == ["GGG", "A", "TTTT"]

    concat = batch.concat()
    assert str(concat) == "ACCGGGTTTT"

    with pytest.raises(ValueError):
        DNABatch([]).concat()


def test_dna_batch_filter_by_len():
    batch = DNABatch([DNA("A"), DNA("CC"), DNA("GGG"), DNA("TTTT")])

    filtered = batch.filter_by_len(2, 3)
    assert [str(s) for s in filtered.to_list()] == ["CC", "GGG"]

    out = batch.filter_by_len(2, 3, inplace=True)
    assert out is None
    assert [str(s) for s in batch.to_list()] == ["CC", "GGG"]

    empty = batch.filter_by_len(5, 2)
    assert len(empty) == 0


def test_dna_batch_count_contains():
    batch = DNABatch([DNA("AAC"), DNA("TT")])

    assert batch.count(DNA("A")) == [2, 0]
    assert batch.contains(DNA("TT")) == [False, True]

    with pytest.raises(TypeError):
        batch.count("A")

    with pytest.raises(TypeError):
        batch.contains("TT")


def test_dna_batch_take_and_slice_errors():
    batch = DNABatch([DNA("A")])

    with pytest.raises(IndexError):
        batch.take([2])

    with pytest.raises(ValueError):
        batch.slice(step=0)


def test_dna_batch_complement_transcribe_translate():
    batch = DNABatch([DNA("ATG"), DNA("CCG")])

    complement = batch.complement()
    assert [str(s) for s in complement.to_list()] == ["TAC", "GGC"]

    out = batch.complement(inplace=True)
    assert out is None
    assert [str(s) for s in batch.to_list()] == ["TAC", "GGC"]

    rna = batch.transcribe()
    assert isinstance(rna, RNABatch)
    assert [str(s) for s in rna.to_list()] == ["UAC", "GGC"]

    protein = DNABatch([DNA("ATG")]).translate()
    assert isinstance(protein, ProteinBatch)
    assert [str(s) for s in protein.to_list()] == ["M"]
