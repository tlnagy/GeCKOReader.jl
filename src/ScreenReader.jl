module ScreenReader

using Distances
using DataStructures
using DataFrames
using Libz
using Bio.Seq
using Formatting
using Gadfly

include("utils.jl")
include("types.jl")
include("io.jl")

export read_seq_file, gen_barcode_mapping, ScreenSample, plot_hit_qualities,
       plot_all_qualities

end # module
