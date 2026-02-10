from logging import raiseExceptions
import pytest

from biorust import DNA, DNABatch, Protein


def test_translate_basic():
    seq = DNA("ATGGCC")
    prot = seq.translate()
    assert isinstance(prot, Protein)
    assert str(prot) == "MA"


def test_translate_lowercase_stop_and_ambiguous():
    assert str(DNA("atgtaa").translate()) == "M*"
    assert str(DNA("ATGNNN").translate()) == "MX"

    # not factor of 3 errors
    with pytest.raises(ValueError):
        DNA("ATGA").translate()


def test_protein_construction():
    p = Protein("M*X")
    assert str(p) == "M*X"

    with pytest.raises(Exception):
        Protein("M#")


# --- Reading frame tests ---


def test_translate_frame_1_drops_trailing():
    seq = DNA("ATGGCCA")
    prot = seq.translate(frame=1)
    assert str(prot) == "MA"


def test_translate_frame_2():
    seq = DNA("AATGGCC")
    prot = seq.translate(frame=2)
    assert str(prot) == "MA"


def test_translate_frame_3():
    seq = DNA("CCATGGCC")
    prot = seq.translate(frame=3)
    assert str(prot) == "MA"


def test_translate_frame_auto():
    # Frame 2 has ATG -> best ORF
    seq = DNA("CATGAAATTT")
    prot = seq.translate(frame="auto")
    assert str(prot) == "MKF"


def test_translate_no_frame_still_strict():
    with pytest.raises(ValueError):
        DNA("ATGA").translate()


def test_translate_frame_invalid():
    seq = DNA("ATG")
    with pytest.raises(ValueError):
        seq.translate(frame=0)
    with pytest.raises(ValueError):
        seq.translate(frame=4)
    with pytest.raises(ValueError):
        seq.translate(frame="bad")


def test_translate_batch_with_frame():
    batch = DNABatch([DNA("ATGGCCA"), DNA("AATGGCC")])
    proteins = batch.translate(frame=1)
    assert proteins.lengths() == [2, 2]
