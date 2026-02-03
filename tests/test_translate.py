import pytest

from biorust import DNA, Protein


def test_translate_basic():
    seq = DNA("ATGGCC")
    prot = seq.translate()
    assert isinstance(prot, Protein)
    assert str(prot) == "MA"


def test_translate_lowercase_stop_and_ambiguous():
    assert str(DNA("atgtaa").translate()) == "M*"
    assert str(DNA("ATGNNN").translate()) == "MX"

    # trailing bases are ignored
    assert str(DNA("ATGA").translate()) == "M"


def test_protein_construction():
    p = Protein("M*X")
    assert str(p) == "M*X"

    with pytest.raises(Exception):
        Protein("M#")
