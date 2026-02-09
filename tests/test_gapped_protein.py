import pytest
from biorust import GappedProtein, Protein


def test_construction_from_str():
    seq = GappedProtein("AC-DE.F")
    assert str(seq) == "AC-DE.F"
    assert len(seq) == 7


def test_invalid_char_rejected():
    with pytest.raises(ValueError, match="invalid character"):
        GappedProtein("AC#DE")


def test_equality():
    a = GappedProtein("AC-DE")
    b = GappedProtein("AC-DE")
    c = GappedProtein("AC-DA")
    assert a == b
    assert a != c


def test_ungapped():
    seq = GappedProtein("A-C.D-E")
    ungapped = seq.ungapped()
    assert isinstance(ungapped, Protein)
    assert str(ungapped) == "ACDE"


def test_slice_returns_gapped_protein():
    seq = GappedProtein("AC-DE")
    sliced = seq[1:4]
    assert isinstance(sliced, GappedProtein)
    assert str(sliced) == "C-D"


def test_index_returns_gapped_protein():
    seq = GappedProtein("AC-DE")
    ch = seq[2]
    assert isinstance(ch, GappedProtein)
    assert str(ch) == "-"


def test_repr():
    seq = GappedProtein("AC-DE")
    assert repr(seq) == 'GappedProtein("AC-DE")'


def test_empty():
    seq = GappedProtein("")
    assert len(seq) == 0
    assert str(seq) == ""


def test_all_iupac_plus_gaps():
    s = "ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz*-."
    seq = GappedProtein(s)
    assert str(seq) == s


def test_negative_index():
    seq = GappedProtein("AC-DE")
    assert str(seq[-1]) == "E"
    assert str(seq[-2]) == "D"
    assert str(seq[-5]) == "A"


def test_negative_slice():
    seq = GappedProtein("AC-DE")
    assert str(seq[-3:]) == "-DE"
    assert str(seq[-4:-1]) == "C-D"


def test_index_out_of_range():
    seq = GappedProtein("AC-DE")
    with pytest.raises(IndexError):
        seq[10]
    with pytest.raises(IndexError):
        seq[-6]
