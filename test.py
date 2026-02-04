from biorust import DNA, complement, DNABatch, DNARecordBatch, read_fasta, read_csv

seq1 = DNA("ATCG")
seq2 = DNA(seq="ATCG")
seq3 = DNA(seq="ATC")
print(seq1 + seq2)
print(seq1 == seq2)
print(seq1 == seq3)
print(seq1 > seq2)
print(seq1[1])
print(seq1[1:3])
print(seq1[-1])
print(seq1 * 3)
print(3 * seq1)

seq2 *= 2
seq3 += seq2
print(2 * seq3)
print(seq1)
print(seq2)

seq1 = DNA("AAAAAAA")
print(seq1.count("AAA"))
print(seq1.count_overlap("AAA"))
seq1 = DNA("ATCTGCATTACG")
print(complement(seq1))
print(seq1.complement())

print(seq1.translate())

batch = DNABatch([DNA("ATGCCGCGTTACGTACG"), DNA("AACGATCGACTACGA")])
batch.reverse_complements(inplace=True)
seq = batch[0]
print(batch[0])
print(batch[0] == seq)

seqs = read_fasta("/Users/nick/Desktop/fasta/sequences.fasta")

seqs = read_csv(
    path="/Users/nick/Library/CloudStorage/Box-Box/Zhang/projects/tf_cloning/Edited-FinalTF sheet.csv",
    id_col="DNASU_Clone_ID",
    seq_col="Sequence",
    on_error="skip"
)

print(len(seqs))
