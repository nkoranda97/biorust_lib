import pytest

from biorust import DNABatch, ProteinBatch, RNA, RNABatch


def test_rna_batch_complement_back_transcribe_translate():
    batch = RNABatch([RNA("AUG"), RNA("CCG")])

    complement = batch.complement()
    assert [str(s) for s in complement.to_list()] == ["UAC", "GGC"]

    out = batch.complement(inplace=True)
    assert out is None
    assert [str(s) for s in batch.to_list()] == ["UAC", "GGC"]

    dna = batch.back_transcribe()
    assert isinstance(dna, DNABatch)
    assert [str(s) for s in dna.to_list()] == ["TAC", "GGC"]

    protein = RNABatch([RNA("AUG")]).translate()
    assert isinstance(protein, ProteinBatch)
    assert [str(s) for s in protein.to_list()] == ["M"]


def test_rna_batch_complement_errors():
    batch = RNABatch([RNA("A")])
    with pytest.raises(TypeError):
        batch.complement(inplace="yes")
