# biorust

For fun biopython and rust-bio inspired library aimed to create a fast bioinformatics library that is a little more explicit and type safe than biopython. Also some SIMD-accelerated alignment, and parallel batch operations. Lots of experimenting with LLM's.

## Install

Install from source

```bash
pip install -e .
```

## Quickstart

```python
from biorust import DNA, Protein, read_fasta

seq = DNA("ATGGCCATTGAATGA")
seq.translate()          # Protein("MAIE*")
seq.reverse_complement() # DNA("TCATTCAATGGCCAT")
seq.gc_content()         # 0.4

records = read_fasta("sequences.fasta")
proteins = records.translate(frame=1)
```

## Sequences

Three core sequence types with full string-like behavior: indexing, slicing, concatenation, comparison, iteration.

```python
from biorust import DNA, RNA, Protein

# DNA
dna = DNA("ATGGCCATTGAATGA")
dna.complement()              # DNA("TACCGGTAACTTACT")
dna.reverse_complement()      # DNA("TCATTCAATGGCCAT")
dna.gc_content()              # 0.4
dna.transcribe()              # RNA("AUGGCCAUUGAAUGA")

# RNA
rna = RNA("AUGGCCAUUGAAUGA")
rna.back_transcribe()         # DNA("ATGGCCATTGAATGA")

# Protein
prot = Protein("MAIE")
prot.molecular_weight()       # MW in Daltons
prot.isoelectric_point()      # predicted pI
prot.hydrophobicity()         # Kyte-Doolittle mean
prot.net_charge(7.4)          # charge at physiological pH
prot.shannon_entropy()        # compositional complexity
```

### String operations

All sequence types support Python string methods like biopython:

```python
dna = DNA("AATGGCCAA")

dna.count("A")                # 4
dna.find("GCC")               # 3
"GCC" in dna                  # True
dna.startswith("AAT")         # True

dna.upper()                   # case conversion
dna.split("GCC")              # [DNA("AAT"), DNA("AA")]
dna.strip("A")                # DNA("TGGCC")

dna[2:5]                      # DNA("TGG") - slicing
dna + DNA("TTT")              # DNA("AATGGCCAATTT") - concatenation
dna * 2                       # DNA("AATGGCCAAAATGGCCAA") - repeat

# However, you can't do this:
dna + "TTT"
```

## Translation

Standard codon table with reading frame support and auto-detection:

```python
dna = DNA("ATGGCCATTGAATGA")
dna.translate()               

# Reading frames - silently drops trailing bases
DNA("ATGGCCA").translate(frame=1)     
DNA("AATGGCC").translate(frame=2)     
DNA("CCATGGCC").translate(frame=3)    

# Auto-detect: picks frame with longest ORF
DNA("CATGAAATTT").translate(frame="auto")  # Protein("MKF")
```

## Primary Sequence Properties

```python
prot = Protein("MKWVTFISLLFLFSSAYS")

prot.molecular_weight()                
prot.isoelectric_point()               
prot.net_charge(7.4)                   
prot.hydrophobicity()                  
prot.hydrophobicity_profile(window=5)  
prot.shannon_entropy()                 

prot.counts()                          
prot.aa_counts_20()                    
prot.aa_frequencies_20()               

prot.has_ambiguous()                   # True if X, B, Z present
prot.unknown_positions()               # indices of ambiguous residues
```

## Batch Operations

Work with sequences in batch:

```python
from biorust import DNABatch, DNA

batch = DNABatch([DNA("ATGGCC"), DNA("ATGAAA"), DNA("ATGTTT")])

batch.lengths()                # [6, 6, 6]
batch.translate()              # ProteinBatch with all translations
batch.translate(frame=1)       # reading frame support on batches too
batch.reverse_complements()    # batch reverse complement
batch.complement()             # batch complement
batch.transcribe()             # DNABatch -> RNABatch

# Filtering and subsetting
batch.filter_by_len(min_len=3, max_len=100)
batch.take([0, 2])            # subset by indices
batch.slice(0, 2)             # subset by range
batch.concat()                # concatenate all into one sequence

# In-place mutations
batch.append(DNA("ATGCCC"))
batch.extend([DNA("GGG")])
batch.reverse_complements(inplace=True)
```

## File I/O

### FASTA

```python
from biorust import read_fasta, write_fasta

# Read - auto-detects DNA/RNA/protein
records = read_fasta("sequences.fasta")
records = read_fasta("sequences.fasta", alphabet="protein")

records.ids()                  # ["seq1", "seq2", ...]
records.descriptions()         # [None, "some desc", ...]
records.seqs()                 # DNABatch of sequences
records[0].seq                 # first sequence
records[0].id                  # first ID

# Write
write_fasta("output.fasta", records)
write_fasta("output.fasta", records, line_width=80)
```

### FASTQ

```python
from biorust import read_fastq, write_fastq

# Read - auto-detects DNA/RNA/protein
records = read_fastq("reads.fastq")
records = read_fastq("reads.fastq", alphabet="dna")

# Write
# Uses a constant quality character because records currently do not store
# per-base quality vectors.
write_fastq("output.fastq", records)
write_fastq("output.fastq", records, quality_char="J")
```

### CSV

Read sequences from a csv by explicitly defining columns

```python
from biorust import read_csv, csv_columns

csv_columns("data.csv")       # ["id", "sequence", "description"]

records = read_csv(
    "data.csv",
    id_col="id",
    seq_col="sequence",
    desc_col="description",
    alphabet="dna",
)

# Skip invalid sequences instead of raising
records = read_csv("data.csv", id_col=0, seq_col=1, on_error="skip")
records.skipped               # [SkippedRecord(row=3, message="..."), ...]
```

## Sequence Records

Records bundle a sequence with metadata:

```python
from biorust import DNARecord, DNARecordBatch, SeqFeature, FeatureLocation

record = DNARecord("gene1", DNA("ATGGCCGAA"), desc="example gene")
record.id                     # "gene1"
record.seq                    # DNA("ATGGCCGAA")
record.description            # "example gene"

# Annotations and features
loc = FeatureLocation(0, 9, strand=1)
feature = SeqFeature("CDS", loc, qualifiers={"gene": ["geneA"]})
record = DNARecord("gene1", DNA("ATGGCCGAA"), features=[feature])

# Record batches
batch = DNARecordBatch([record1, record2, record3])
batch.translate()              # ProteinRecordBatch (preserves IDs)
batch.filter_empty()           # remove empty sequences
batch.reverse_complements()    # features update coordinates automatically
```

## Pairwise Alignment

Smith-Waterman (local) and Needleman-Wunsch (global) with SIMD acceleration:

```python
from biorust import DNA, Protein, Scoring, align_local, align_global

# Simple scoring
scoring = Scoring(match_score=2, mismatch_score=-1, gap_open=-5, gap_extend=-1)

# Or use a substitution matrix
scoring = Scoring(matrix="BLOSUM62", gap_open=-10, gap_extend=-1)

# Score-only (SIMD-accelerated)
result = align_local(DNA("ACGTACGT"), DNA("ACGTAAACGT"), scoring)
result.score                   # alignment score

# With traceback
result = align_global(
    DNA("ACGTACGT"),
    DNA("ACGTAAACGT"),
    scoring,
    traceback=True,
)
result.score
result.cigar                   # [("M", 4), ("I", 2), ("M", 4)]
result.aligned_strings()       # ("ACGT--ACGT", "ACGTAAACGT")
print(result.alignment_diagram())

# Protein alignment
result = align_local(
    Protein("MVLSPADKTNVK"),
    Protein("MVHLTPEEKSAVT"),
    Scoring(matrix="BLOSUM62", gap_open=-10, gap_extend=-1),
    traceback=True,
)
```

## Multiple Sequence Alignment

Via Clustal Omega (requires `clustalo` binary):

```python
from biorust import read_fasta, msa

records = read_fasta("sequences.fasta")
alignment = msa(records, algorithm="clustalo")

alignment.width                # alignment length
len(alignment)                 # number of sequences
alignment.ids()                # sequence IDs

print(alignment)               # Clustal-style diagram with conservation
```

## Phylogenetics

Build trees from multiple sequence alignments:

```python
from biorust import distance_matrix, build_tree

# Distance matrix from alignment
dm = distance_matrix(alignment, model="jc69")  # DNA: p-distance, jc69, k2p
dm = distance_matrix(alignment, model="poisson")  # Protein: p-distance, poisson

dm.labels()                    # sequence labels
dm.get(0, 1)                  # distance between sequences 0 and 1
dm.to_list_of_lists()         # full matrix as 2D list

# Build tree
tree = build_tree(dm, method="nj")      # neighbor-joining (unrooted)
tree = build_tree(dm, method="upgma")   # UPGMA (rooted)

tree.to_newick()               # Newick format string
tree.leaf_labels()             # taxon names
tree.is_rooted()               # True for UPGMA, False for NJ
print(tree.ascii_diagram())    # ASCII tree visualization
```
