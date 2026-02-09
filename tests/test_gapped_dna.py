import pytest
from biorust import GappedDNA, DNA


def test_construction_from_str():
    seq = GappedDNA("AC-GT.N")
    assert str(seq) == "AC-GT.N"
    assert len(seq) == 7


def test_invalid_char_rejected():
    with pytest.raises(ValueError, match="invalid character"):
        GappedDNA("AC#GT")


def test_equality():
    a = GappedDNA("AC-GT")
    b = GappedDNA("AC-GT")
    c = GappedDNA("AC-GA")
    assert a == b
    assert a != c


def test_ungapped():
    seq = GappedDNA("A-C.G-T")
    ungapped = seq.ungapped()
    assert isinstance(ungapped, DNA)
    assert str(ungapped) == "ACGT"


def test_slice_returns_gapped_dna():
    seq = GappedDNA("AC-GT")
    sliced = seq[1:4]
    assert isinstance(sliced, GappedDNA)
    assert str(sliced) == "C-G"


def test_index_returns_gapped_dna():
    seq = GappedDNA("AC-GT")
    ch = seq[2]
    assert isinstance(ch, GappedDNA)
    assert str(ch) == "-"


def test_repr():
    seq = GappedDNA("AC-GT")
    assert repr(seq) == 'GappedDNA("AC-GT")'


def test_empty():
    seq = GappedDNA("")
    assert len(seq) == 0
    assert str(seq) == ""


def test_all_iupac_plus_gaps():
    s = "ACGTRYSWKMBDHVNacgtryswkmbdhvn-."
    seq = GappedDNA(s)
    assert str(seq) == s


def test_negative_index():
    seq = GappedDNA("AC-GT")
    assert str(seq[-1]) == "T"
    assert str(seq[-2]) == "G"
    assert str(seq[-5]) == "A"


def test_negative_slice():
    seq = GappedDNA("AC-GT")
    assert str(seq[-3:]) == "-GT"
    assert str(seq[-4:-1]) == "C-G"


def test_index_out_of_range():
    seq = GappedDNA("AC-GT")
    with pytest.raises(IndexError):
        seq[10]
    with pytest.raises(IndexError):
        seq[-6]
