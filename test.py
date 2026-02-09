"""biorust showcase script.

Run:
  uv run python test.py

This script demonstrates:
- DNA/RNA/Protein basic operations
- Features + annotations on records
- Batch operations and reverse complements
- FASTA read/write (local temp files)
- CSV read with error handling
- Global/local alignments with traceback
"""

from __future__ import annotations

from pathlib import Path
import tempfile

from biorust import (
    AlignmentResult,
    DNA,
    DNABatch,
    DNARecord,
    DNARecordBatch,
    FeatureLocation,
    RNABatch,
    RNA,
    RNARecord,
    Scoring,
    SeqFeature,
    align_global,
    align_local,
    complement,
    read_csv,
    read_fasta,
    write_fasta,
    msa_clustalo
)


def header(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def demo_sequences() -> None:
    header("1) Sequence Basics")

    dna = DNA("ATCG")
    dna2 = DNA("ATCG")
    dna3 = DNA("ATC")

    print("DNA:", dna)
    print("Equality:", dna == dna2, dna == dna3)
    print("Ordering:", dna > dna2)
    print("Indexing:", dna[1], dna[-1])
    print("Slicing:", dna[1:3])
    print("Repeat:", dna * 3)
    print("Count:", DNA("AAAAAAA").count("AAA"))
    print("Count overlap:", DNA("AAAAAAA").count_overlap("AAA"))
    print("Complement (function):", complement(dna))
    print("Complement (method):", dna.complement())

    rna = dna.transcribe()
    print("Transcribe DNA -> RNA:", rna)
    print("Back-transcribe RNA -> DNA:", rna.back_transcribe())

    protein = DNA("ATGGCC").translate()
    print("Translate DNA -> Protein:", protein)
    print("Protein length:", len(protein))


def demo_features_annotations() -> None:
    header("2) Features + Annotations")

    loc = FeatureLocation(start=1, end=4, strand=1)
    feat = SeqFeature("gene", loc, {"gene": ["abc"]})

    rec = DNARecord(
        "seq1",
        DNA("ATGCAT"),
        desc="example record",
        features=[feat],
        annotations={"source": ["showcase"], "note": ["demo"]},
    )

    rna_rec = RNARecord(
        "rna1",
        RNA("AUGCUU"),
        features=[feat],
        annotations={"source": ["showcase"]},
    )

    print("Record:", rec)
    print("RNARecord:", rna_rec)
    print("Feature type:", rec.features[0].feature_type)
    print(
        "Feature location:",
        rec.features[0].location.start,
        rec.features[0].location.end,
        rec.features[0].location.strand,
    )
    print("Qualifiers:", rec.features[0].qualifiers)
    print("Annotations:", rec.annotations)


def demo_batches() -> None:
    header("3) Batch Operations")

    batch = DNABatch([DNA("ATGCCG"), DNA("AACGAT")])
    print("Lengths:", batch.lengths())
    print("Reverse complements (in-place)")
    batch.reverse_complements(inplace=True)
    print("Batch[0]:", batch[0])

    rna_batch = RNABatch([RNA("AUGC"), RNA("UUAA")])
    print("RNA batch lengths:", rna_batch.lengths())


def demo_fasta_io(tmp_dir: Path) -> None:
    header("4) FASTA Read/Write")

    records = [
        DNARecord("seq1", DNA("ACGT"), "some desc"),
        DNARecord("seq2", DNA("A")),
    ]

    fasta_path = tmp_dir / "demo.fasta"
    write_fasta(str(fasta_path), records, line_width=2)
    print("Wrote:", fasta_path)

    batch = read_fasta(str(fasta_path))
    print("Read batch size:", len(batch))
    print("First id/desc/seq:", batch[0].id, batch[0].description, batch[0].seq)


def demo_csv_io(tmp_dir: Path) -> None:
    header("5) CSV Read")

    csv_path = tmp_dir / "demo.csv"
    csv_path.write_text("id,Sequence,Note\ns1,ATGC,ok\ns2,ATG#BAD,skip\ns3,AACG,ok\n")

    batch = read_csv(
        path=str(csv_path),
        id_col="id",
        seq_col="Sequence",
        desc_col="Note",
        on_error="skip",
    )

    print("CSV batch size:", len(batch))
    print("Skipped:", [s.row for s in batch.skipped])
    print("First record:", batch[0].id, batch[0].description, batch[0].seq)


def demo_alignment() -> None:
    header("6) Alignment")

    seq1 = DNA("ACGTACGT")
    seq2 = DNA("ACGTTACGT")
    scoring = Scoring(gap_open=-10, gap_extend=-0.5)

    global_res: AlignmentResult = align_global(
        seq1, seq2, scoring=scoring, traceback=True
    )
    local_res: AlignmentResult = align_local(
        seq1, seq2, scoring=scoring, traceback=True
    )

    print("Global score:", global_res.score)
    print("Local score:", local_res.score)

    aligned1, aligned2 = global_res.aligned_strings()
    print("\nGlobal aligned sequences:")
    print(aligned1)
    print(aligned2)
    print("\nAlignment diagram:\n" + global_res.alignment_diagram())


def main() -> None:
    demo_sequences()
    demo_features_annotations()
    demo_batches()

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        demo_fasta_io(tmp_dir)
        demo_csv_io(tmp_dir)

    demo_alignment()

    records = DNARecordBatch([
        DNARecord("seq1", DNA("ATGCTAGCTAG")),
        DNARecord("seq2", DNA("ATGCGAGCTAG")),
        DNARecord("seq3", DNA("ATGCTAGATAG")),
    ])

    msa = msa_clustalo(records)
    print("\n")
    print(msa)

    


if __name__ == "__main__":
    main()
