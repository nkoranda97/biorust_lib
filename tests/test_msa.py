import pytest
from biorust import AlignmentDNA, AlignmentProtein, GappedDNA, GappedProtein


def _make_alignment():
    return AlignmentDNA(
        [
            ("seq1", GappedDNA("AC-GT")),
            ("seq2", GappedDNA("ACAGT")),
            ("seq3", GappedDNA("AC-GA")),
        ]
    )


def test_rejects_unequal_lengths():
    with pytest.raises(ValueError, match="equal length"):
        AlignmentDNA(
            [
                ("a", GappedDNA("ACGT")),
                ("b", GappedDNA("AC")),
            ]
        )


def test_rejects_empty():
    with pytest.raises(ValueError, match="at least one"):
        AlignmentDNA([])


def test_len_and_width():
    aln = _make_alignment()
    assert len(aln) == 3
    assert aln.width == 5


def test_ids_preserved():
    aln = _make_alignment()
    assert aln.ids() == ["seq1", "seq2", "seq3"]


def test_seqs_returns_gapped_dna():
    aln = _make_alignment()
    seqs = aln.seqs()
    assert len(seqs) == 3
    assert all(isinstance(s, GappedDNA) for s in seqs)
    assert str(seqs[0]) == "AC-GT"


def test_getitem_int():
    aln = _make_alignment()
    id_, seq = aln[0]
    assert id_ == "seq1"
    assert isinstance(seq, GappedDNA)
    assert str(seq) == "AC-GT"


def test_getitem_negative():
    aln = _make_alignment()
    id_, seq = aln[-1]
    assert id_ == "seq3"
    assert str(seq) == "AC-GA"


def test_getitem_out_of_range():
    aln = _make_alignment()
    with pytest.raises(IndexError):
        aln[10]


def test_getitem_slice():
    aln = _make_alignment()
    sub = aln[0:2]
    assert isinstance(sub, AlignmentDNA)
    assert len(sub) == 2
    assert sub.ids() == ["seq1", "seq2"]
    assert sub.width == 5


def test_getitem_slice_negative():
    aln = _make_alignment()
    sub = aln[-2:]
    assert len(sub) == 2
    assert sub.ids() == ["seq2", "seq3"]


def test_getitem_slice_step():
    aln = _make_alignment()
    sub = aln[::2]
    assert len(sub) == 2
    assert sub.ids() == ["seq1", "seq3"]


def test_repr():
    aln = _make_alignment()
    assert repr(aln) == "AlignmentDNA(n=3, width=5)"


# ---- AlignmentProtein tests ----


def _make_protein_alignment():
    return AlignmentProtein(
        [
            ("seq1", GappedProtein("AC-DE")),
            ("seq2", GappedProtein("ACADE")),
            ("seq3", GappedProtein("AC-DA")),
        ]
    )


def test_protein_rejects_unequal_lengths():
    with pytest.raises(ValueError, match="equal length"):
        AlignmentProtein(
            [
                ("a", GappedProtein("ACDE")),
                ("b", GappedProtein("AC")),
            ]
        )


def test_protein_rejects_empty():
    with pytest.raises(ValueError, match="at least one"):
        AlignmentProtein([])


def test_protein_len_and_width():
    aln = _make_protein_alignment()
    assert len(aln) == 3
    assert aln.width == 5


def test_protein_ids_preserved():
    aln = _make_protein_alignment()
    assert aln.ids() == ["seq1", "seq2", "seq3"]


def test_protein_seqs_returns_gapped_protein():
    aln = _make_protein_alignment()
    seqs = aln.seqs()
    assert len(seqs) == 3
    assert all(isinstance(s, GappedProtein) for s in seqs)
    assert str(seqs[0]) == "AC-DE"


def test_protein_getitem_int():
    aln = _make_protein_alignment()
    id_, seq = aln[0]
    assert id_ == "seq1"
    assert isinstance(seq, GappedProtein)
    assert str(seq) == "AC-DE"


def test_protein_getitem_negative():
    aln = _make_protein_alignment()
    id_, seq = aln[-1]
    assert id_ == "seq3"
    assert str(seq) == "AC-DA"


def test_protein_getitem_out_of_range():
    aln = _make_protein_alignment()
    with pytest.raises(IndexError):
        aln[10]


def test_protein_getitem_slice():
    aln = _make_protein_alignment()
    sub = aln[0:2]
    assert isinstance(sub, AlignmentProtein)
    assert len(sub) == 2
    assert sub.ids() == ["seq1", "seq2"]
    assert sub.width == 5


def test_protein_repr():
    aln = _make_protein_alignment()
    assert repr(aln) == "AlignmentProtein(n=3, width=5)"
