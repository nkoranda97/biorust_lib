from biorust import DNA, DNABatch, DNARecord, DNARecordBatch


def test_dna_record_batch_basic():
    r1 = DNARecord("id1", DNA("ATGC"), "desc1")
    r2 = DNARecord("id2", DNA("AACG"), None)

    batch = DNARecordBatch([r1, r2])
    assert len(batch) == 2
    assert batch.ids() == ["id1", "id2"]
    assert batch.descriptions() == ["desc1", None]
    assert isinstance(batch.seqs(), DNABatch)
    assert [str(s) for s in batch.seqs().to_list()] == ["ATGC", "AACG"]


def test_dna_record_batch_slice_alignment():
    r1 = DNARecord("id1", DNA("ATGC"), "desc1")
    r2 = DNARecord("id2", DNA("AACG"), None)
    r3 = DNARecord("id3", DNA("TTAA"), "desc3")
    batch = DNARecordBatch([r1, r2, r3])

    sub = batch[1:3]
    assert sub.ids() == ["id2", "id3"]
    assert sub.descriptions() == [None, "desc3"]
    assert [str(s) for s in sub.seqs().to_list()] == ["AACG", "TTAA"]


def test_dna_record_batch_reverse_complements():
    r1 = DNARecord("id1", DNA("ATGC"), "desc1")
    r2 = DNARecord("id2", DNA("AACG"), None)
    batch = DNARecordBatch([r1, r2])

    out = batch.reverse_complements()
    assert out is not batch
    assert out.ids() == ["id1", "id2"]
    assert out.descriptions() == ["desc1", None]
    assert [str(s) for s in out.seqs().to_list()] == ["GCAT", "CGTT"]

    inplace = batch.reverse_complements(inplace=True)
    assert inplace is None
    assert batch.ids() == ["id1", "id2"]
    assert batch.descriptions() == ["desc1", None]
    assert [str(s) for s in batch.seqs().to_list()] == ["GCAT", "CGTT"]
