immutable ScreenSample
    "name of original file"
    filename::String
    "total number of reads"
    reads::Int
    "Number of called hits"
    hits::Int
    "Number of high quality reads that didn't match"
    no_match::Int
    "Number of reads with too many errors in their constant region"
    constant_error::Int
    "Number of reads with low quality sgRNA regions"
    low_quality::Int
    "Number of reads with a priming site starting at each position of the read"
    alt_staggers::Vector{Int}
    counts::Dict{String, Int}
    hit_qualities::Matrix{Int}
    all_qualities::Matrix{Int}

    function ScreenSample(filename, read_length, max_phred_quality)
        new(filename, 0, 0, 0, 0, 0, zeros(Int, read_length),
            Dict{String, Int}(),
            zeros(Int, max_phred_quality, read_length),
            zeros(Int, max_phred_quality, read_length))
    end
end


function Base.show(io::IO, s::ScreenSample)
    println(io, typeof(s))
    println(io, "\tFilename: $(s.filename)")
    println(io, "\tHits: ", format(s.hits, commas=true))
    println(io, "\tReads: ", format(s.reads, commas=true))
    println(io, "\t% reads matched: $(round(s.hits/s.reads*100, 4))%")
    println(io, "\t% unmatchable reads: $(round(s.no_match/(s.hits+s.no_match)*100, 4))%")
    println(io, "\t% reads had errors in the constant region: $(round(s.constant_error/s.reads*100, 4))%")
    println(io, "\t% reads had low quality sgrna regions: $(round(s.low_quality/s.reads*100, 4))%")
end
