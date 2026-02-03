from biorust import DNA

seq1 = DNA("ATCG")
seq2 = DNA(seq="ATCG")
seq3 = DNA(seq="ATC")
print(seq1 + "ATTT")
print(seq1 + seq2)
print(seq1 == seq2)
print(seq1 == seq3)
print(seq1 > seq2)
print(seq1[1])
print(seq1[1:3])
print(seq1[-1])