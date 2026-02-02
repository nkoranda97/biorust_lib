from biorust import DNA

seq = DNA("ATCG")
print(seq + "ATTT")
print(seq + DNA("ATTT"))