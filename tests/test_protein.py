import pytest

from biorust import Protein


def test_protein_basic_ops():
    seq1 = Protein("ACD")
    seq2 = Protein("ACD")
    seq3 = Protein("AC")

    # addition
    assert str(seq1 + "EF") == "ACDEF"
    assert str(seq1 + seq2) == "ACDACD"

    # equality + ordering
    assert (seq1 == seq2) is True
    assert (seq1 == seq3) is False
    assert (seq1 > seq2) is False

    # indexing / slicing return Protein (not bytes)
    one = seq1[1]
    assert isinstance(one, Protein)
    assert str(one) == "C"

    sl = seq1[1:3]
    assert isinstance(sl, Protein)
    assert str(sl) == "CD"

    last = seq1[-1]
    assert isinstance(last, Protein)
    assert str(last) == "D"

    # multiplication
    assert str(seq1 * 2) == "ACDACD"
    assert str(2 * seq1) == "ACDACD"


def test_protein_count_contains_find():
    seq = Protein("AAAAA")

    assert seq.count("AA") == 2
    assert seq.count_overlap("AA") == 4

    assert "A" in seq
    assert "Z" not in seq

    assert seq.find("AAA") == 0
    assert seq.rfind("AAA") == 2


def test_protein_split_strip_case():
    seq = Protein("AAACCCAAA")

    assert [str(p) for p in seq.split("CCC")] == ["AAA", "AAA"]
    assert [str(p) for p in seq.rsplit("A", 1)] == ["AAACCCAA", ""]

    assert str(seq.strip("A")) == "CCC"
    assert str(seq.lstrip("A")) == "CCCAAA"
    assert str(seq.rstrip("A")) == "AAACCC"

    mixed = Protein("aCdE")
    assert str(mixed.upper()) == "ACDE"
    assert str(mixed.lower()) == "acde"


def test_protein_startswith_endswith():
    seq = Protein("ACDE")

    assert seq.startswith("AC")
    assert seq.startswith(("ZZ", "AC"))
    assert seq.startswith(65)  # "A"

    assert seq.endswith("DE")
    assert seq.endswith(("ZZ", "DE"))
    assert seq.endswith(69)  # "E"

    with pytest.raises(Exception):
        Protein("A#")
