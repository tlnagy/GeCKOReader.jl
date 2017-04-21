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
