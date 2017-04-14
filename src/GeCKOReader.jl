module GeCKOReader

using Distances
using DataStructures
using DataFrames
using Libz
using Bio.Seq
using Formatting
using Gadfly
using Query
using HypothesisTests

include("utils.jl")
include("types.jl")
include("io.jl")
include("processing.jl")

export read_seq_file, gen_barcode_mapping, ScreenSample, plot_hit_qualities,
       plot_all_qualities, get_rel_freqs, get_log2fc_pvals

end # module
