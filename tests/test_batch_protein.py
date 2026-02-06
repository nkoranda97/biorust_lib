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


def test_protein_batch_biology_metrics():
    batch = ProteinBatch([Protein("ACD")])

    reversed_batch = batch.reverse()
    assert [str(s) for s in reversed_batch.to_list()] == ["DCA"]

    out = batch.reverse(inplace=True)
    assert out is None
    assert [str(s) for s in batch.to_list()] == ["DCA"]

    counts = batch.counts()
    assert counts == [[("A", 1), ("C", 1), ("D", 1)]]

    frequencies = batch.frequencies()
    assert len(frequencies) == 1
    assert abs(sum(val for _, val in frequencies[0]) - 1.0) < 1e-9

    aa_counts_20 = batch.aa_counts_20()
    assert len(aa_counts_20) == 1
    assert len(aa_counts_20[0]) == 20
    assert aa_counts_20[0][0] == ("A", 1)

    aa_frequencies_20 = batch.aa_frequencies_20()
    assert len(aa_frequencies_20) == 1
    assert len(aa_frequencies_20[0]) == 20

    entropy = batch.shannon_entropy()
    assert len(entropy) == 1
    assert entropy[0] > 0.0

    weights = batch.molecular_weight()
    assert len(weights) == 1
    assert weights[0] > 0.0

    hydrophobicity = batch.hydrophobicity()
    assert len(hydrophobicity) == 1

    profile = batch.hydrophobicity_profile(2)
    assert len(profile) == 1
    assert len(profile[0]) == 2

    charges = batch.net_charge(7.0)
    assert len(charges) == 1

    pis = batch.isoelectric_point()
    assert len(pis) == 1

    amb = ProteinBatch([Protein("AX"), Protein("AC")])
    assert amb.has_ambiguous() == [True, False]
    assert amb.unknown_positions() == [[1], []]
