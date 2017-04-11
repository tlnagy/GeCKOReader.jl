using GeCKOReader
using Base.Test
using Bio.Seq

# 1 error
seq_a = dna"ATCGT"
seq_b = dna"ATAGT"
@test GeCKOReader.hamming(seq_a, seq_b) == 1

# 0 error
seq_a = dna"ATCGT"
seq_b = dna"ATCGT"
@test GeCKOReader.hamming(seq_a, seq_b) == 0

# All error
seq_a = dna"AAAAA"
seq_b = dna"TTTTT"
@test GeCKOReader.hamming(seq_a, seq_b) == length(seq_a)

# Empty sequence test
seq_a = dna""
seq_b = dna""
@test GeCKOReader.hamming(seq_a, seq_b) == 0

# Comparing DNA sequences of different lengths is undefined
seq_a = dna"CGA"
seq_b = dna"CGAA"
@test_throws DimensionMismatch GeCKOReader.hamming(seq_a, seq_b) == 0
