import pytest

from biorust import DNA


def test_construction_and_basic_ops():
    seq1 = DNA("ATCG")
    seq2 = DNA(seq="ATCG")
    seq3 = DNA(seq="ATC")

    # addition
    assert str(seq1 + "ATTT") == "ATCGATTT"
    assert str(seq1 + seq2) == "ATCGATCG"

    # equality + ordering
    assert (seq1 == seq2) is True
    assert (seq1 == seq3) is False
    assert (seq1 > seq2) is False

    # indexing / slicing return DNA (not bytes)
    one = seq1[1]
    assert isinstance(one, DNA)
    assert str(one) == "T"

    sl = seq1[1:3]
    assert isinstance(sl, DNA)
    assert str(sl) == "TC"

    last = seq1[-1]
    assert isinstance(last, DNA)
    assert str(last) == "G"

    # multiplication
    assert str(seq1 * 3) == "ATCGATCGATCG"
    assert str(3 * seq1) == "ATCGATCGATCG"

    seq2 *= 2
    assert str(seq2) == "ATCGATCG"

    seq3 += seq2
    assert str(seq3) == "ATCATCGATCG"

    assert str(2 * seq3) == "ATCATCGATCGATCATCGATCG"


def test_count_and_count_overlap():
    seq = DNA("AAAAAAA")

    # non-overlapping occurrences: AAA|AAA|A => 2
    assert seq.count("AAA") == 2

    # overlapping occurrences: positions 0..4 => 5
    assert seq.count_overlap("AAA") == 5


def test_invalid_input_rejected():
    with pytest.raises(Exception):
        DNA("ATC#")
