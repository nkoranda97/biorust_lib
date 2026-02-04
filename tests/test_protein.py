import pytest

from biorust import Protein


def test_protein_basic_ops():
    seq1 = Protein("ACD")
    seq2 = Protein("ACD")
    seq3 = Protein("AC")

    # addition
    assert str(seq1 + Protein("EF")) == "ACDEF"
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

    with pytest.raises(TypeError):
        seq1 + "EF"

    with pytest.raises(TypeError):
        "EF" + seq1


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


def test_protein_properties():
    seq = Protein("ACAC")
    assert str(seq.reverse()) == "CACA"

    counts = dict(seq.counts())
    assert counts["A"] == 2
    assert counts["C"] == 2

    aa_counts = dict(seq.aa_counts_20())
    assert aa_counts["A"] == 2
    assert aa_counts["C"] == 2

    freq = dict(seq.aa_frequencies_20())
    assert freq["A"] == pytest.approx(0.5)
    assert freq["C"] == pytest.approx(0.5)

    assert seq.shannon_entropy() == pytest.approx(1.0)

    mw = Protein("AC").molecular_weight()
    assert mw == pytest.approx(192.23288, rel=1e-6)

    hydro = Protein("ACD").hydrophobicity()
    assert hydro == pytest.approx(0.2666666667, rel=1e-6)
    profile = Protein("ACD").hydrophobicity_profile(2)
    assert profile == pytest.approx([2.15, -0.5], rel=1e-6)

    pi = Protein("AC").isoelectric_point()
    charge_pi = Protein("AC").net_charge(pi)
    assert pi >= 0.0 and pi <= 14.0
    assert charge_pi == pytest.approx(0.0, abs=1e-2)


def test_protein_ambiguity_helpers():
    seq = Protein("ACBX")
    assert seq.has_ambiguous() is True
    assert seq.unknown_positions() == [2, 3]

    with pytest.raises(ValueError):
        seq.validate_strict_20()
