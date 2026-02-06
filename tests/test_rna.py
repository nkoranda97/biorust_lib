import pytest

from biorust import DNA, RNA, Protein


def test_rna_basic_ops():
    rna = RNA("AUGC")
    assert str(rna) == "AUGC"
    assert str(rna.complement()) == "UACG"
    assert str(rna.reverse_complement()) == "GCAU"


def test_rna_back_transcribe():
    rna = RNA("AUGC")
    dna = rna.back_transcribe()
    assert isinstance(dna, DNA)
    assert str(dna) == "ATGC"


def test_dna_transcribe():
    dna = DNA("ATGC")
    rna = dna.transcribe()
    assert isinstance(rna, RNA)
    assert str(rna) == "AUGC"


def test_rna_translate():
    rna = RNA("AUGGCC")
    prot = rna.translate()
    assert isinstance(prot, Protein)
    assert str(prot) == "MA"


def test_rna_invalid_char():
    with pytest.raises(Exception):
        RNA("ATGC")
