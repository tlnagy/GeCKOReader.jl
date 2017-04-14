using GeCKOReader
using Base.Test
using Bio.Seq
using Query
using HypothesisTests
using DataFrames

@testset "Hamming" begin

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

end

@testset "Frequencies" begin
    # generate fake sample
    testsample = ScreenSample("test", 51, 41)
    testsample.counts["HGLibA_64614"] = 20
    testsample.counts["HGLibA_64466"] = 30
    testsample.counts["HGLibA_65230"] = 25
    testsample.counts["HGLibA_54150"] = 50

    callable_bc = readtable("callable_sgrnas.csv")

    freqs = GeCKOReader.get_rel_freqs(testsample, callable_bc)

    trimmed_freqs = @from i in freqs begin
        @where i.freqs > 0
        @select i
        @collect DataFrame
    end

    @test sum(trimmed_freqs[:freqs]) == 1

    @test all(isapprox.(Array(trimmed_freqs[:rel_freqs]), [2.0, 1.2, 0.8, 1.0]))
end

@testset "Log2 FC" begin
    df = DataFrame(gene=["A", "B", "A", "A", "B", "A", "B", "C", "C"],
                   log2fc=[5.0, 0.375, 3.0, 6.0, 1.5, Inf, -Inf, Inf, Inf])
    negcontrols = [2.58155,2.29046,1.66488,-0.149982,1.67709]

    result = GeCKOReader.compute_gene_stats(df, negcontrols)

    @test size(result) == (2, 5)

    @test result[1, :med_log2fc] == 5.0
    @test result[2, :med_log2fc] == 0.9375

    pval = -log10(pvalue(MannWhitneyUTest([5.0, 3.0, 6.0], negcontrols)))

    @test result[1, :pval] == pval
end
