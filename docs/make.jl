using Documenter
using GeCKOReader
using Bio

makedocs(
    format = :html,
    sitename = "GeCKOReader.jl",
    modules = [GeCKOReader],
    pages = [
        "index.md"
    ]
)

# deploy off of julia 0.5 instead of nightly
deploydocs(
    julia  = "0.5"
)
