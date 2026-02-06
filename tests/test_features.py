import pytest

from biorust import DNA, DNARecord, DNARecordBatch, FeatureLocation, SeqFeature


def test_feature_location_validation():
    with pytest.raises(Exception):
        FeatureLocation(5, 2)
    with pytest.raises(Exception):
        FeatureLocation(0, 1, 0)


def test_seq_feature_and_record_metadata():
    loc = FeatureLocation(1, 3, 1)
    feat = SeqFeature("gene", loc, {"gene": ["abc"]})
    rec = DNARecord(
        "id1",
        DNA("ATGC"),
        features=[feat],
        annotations={"source": ["test"]},
    )

    assert rec.features[0].feature_type == "gene"
    assert rec.features[0].location.start == 1
    assert rec.features[0].location.end == 3
    assert rec.features[0].location.strand == 1
    assert rec.features[0].qualifiers == {"gene": ["abc"]}
    assert rec.annotations == {"source": ["test"]}

    batch = DNARecordBatch([rec])
    assert batch[0].features[0].feature_type == "gene"
