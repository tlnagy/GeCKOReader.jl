module GeCKOReader

using Distances
using DataStructures
using DataFrames
using Libz
using BioSequences
using Formatting
using Gadfly
using Query
using HypothesisTests

include("utils.jl")
include("types.jl")
include("io.jl")
include("processing.jl")
include("plotting.jl")

export read_seq_file,
       gen_barcode_mapping,
       ScreenSample,
       get_rel_freqs,
       get_log2fc_pvals,
       fit_lines,
       calc_log2fc,
       compute_gene_stats

end # module
