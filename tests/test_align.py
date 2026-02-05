import pytest

from biorust import AlignmentResult, DNA, Protein, Scoring, align_global, align_local


def test_align_local_dna_traceback():
    q = DNA("ACGT")
    t = DNA("ACGT")
    scoring = Scoring(gap_open=-2.0, gap_extend=-1.0)
    res = align_local(q, t, scoring, traceback=True)
    assert isinstance(res, AlignmentResult)
    assert res.score == 20.0
    assert res.cigar is not None
    assert res.cigar == [("M", 4)]
    q_str, mid, t_str = res.aligned_strings(q, t)
    assert q_str == "ACGT"
    assert t_str == "ACGT"
    assert mid == "||||"


def test_align_global_dna_traceback():
    q = DNA("ACGT")
    t = DNA("ACG")
    scoring = Scoring(gap_open=-2.0, gap_extend=-1.0)
    res = align_global(q, t, scoring, traceback=True)
    assert res.score == 13.0
    assert res.cigar is not None
    assert sum(n for _, n in res.cigar) == 4
    q_str, mid, t_str = res.aligned_strings(q, t)
    assert q_str.replace("-", "") == "ACGT"
    assert t_str.replace("-", "") == "ACG"
    assert len(q_str) == len(t_str) == len(mid)


def test_align_type_mismatch():
    q = DNA("ACGT")
    t = Protein("ACGT")
    scoring = Scoring()
    with pytest.raises(ValueError):
        align_local(q, t, scoring, traceback=False)


def test_align_protein():
    q = Protein("ACDE")
    t = Protein("ACDE")
    scoring = Scoring(
        match_score=2,
        mismatch_score=-1,
        gap_open=-3.0,
        gap_extend=-1.0,
        matrix=None,
        use_matrix=False,
    )
    res = align_global(q, t, scoring, traceback=False)
    assert res.score == 8.0
