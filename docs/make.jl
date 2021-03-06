using Documenter
using GeCKOReader
using Bio

makedocs(
    format = :html,
    sitename = "GeCKOReader.jl",
    modules = [GeCKOReader],
    pages = [
        "Home" => "index.md"
    ]
)

# deploy off of julia 0.5 instead of nightly
deploydocs(
    repo="github.com/tlnagy/GeCKOReader.jl.git",
    julia  = "0.5",
    deps = nothing,
    make = nothing,
    target = "build"
)
