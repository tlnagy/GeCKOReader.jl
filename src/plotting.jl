module Plot

using Gadfly
using Query
using DataFrames
using GeCKOReader

export plot_all_qualities, 
       plot_hit_qualities,
       plot_guide_trajectories

function _plot(s::Matrix{Int})
    rows = Int[]
    trim_result = []
    for row in 1:size(s)[1]
        if sum(s[row, :] .!= 0) > 0
            push!(trim_result, s[row, :])
            push!(rows, row)
        end
    end
    trim_result = transpose(hcat(trim_result...))
    spy(trim_result, Scale.y_discrete(labels = i ->rows[i]),
        Guide.xticks(ticks=collect(1:3:size(trim_result)[2])))
end


"""
    plot_all_qualities(s::ScreenSample)

Generate a 2D histogram of read qualities from the FASTQ file per position in
the read using **all** of the read info.
"""
plot_all_qualities(s::ScreenSample) = _plot(s.all_qualities)


"""
    plot_hit_qualities(s::ScreenSample)

Generate a 2D histogram of read qualities from the FASTQ file per position in
the read using **only reads that were matched**.
"""
plot_hit_qualities(s::ScreenSample) = _plot(s.hit_qualities)


"""
    plot_guide_trajectories(lines::DataFrame, gene::String)

Plot the trajectories of the guides for a specific `gene` given the data in
`lines`. The `lines` dataframe should be formatted like the output of
[`fit_lines`](@ref Main.GeCKOReader.fit_lines).
"""
function plot_guide_trajectories(lines::DataFrame, gene::String)

    (!([:t0_rel_freqs, :t1_rel_freqs, :t2_rel_freqs,
        :intercept, :slope, :id, :gene] ⊆ names(lines))) && error("DataFrame malformed missing required columns")

    df = @from i in lines begin
        @where i.gene == gene
        @select i
        @collect DataFrame
    end

    lens = size(df, 1)
    tmp = DataFrame(x=repeat([0,1,2], inner=lens),
                    y=vcat(df[:t0_rel_freqs], df[:t1_rel_freqs], df[:t2_rel_freqs]),
                    guides=repeat(df[:id], outer=3))

    xs = repeat([0, 2], inner=lens)
    tmp2 = DataFrame(x= xs,
                     y=vcat(df[:intercept], df[:intercept] + df[:slope].*2),
                     guides=repeat(df[:id], outer=2))

    plot(
        layer(tmp, x=:x, y=:y, color=:guides, Geom.point),
        layer(tmp2, x=:x, y=:y, color=:guides, Geom.line),
        Guide.xlabel("# of osmoshock treatments"),
        Guide.ylabel("Relative frequencies of guide"),
        Guide.title(gene),
        Guide.xticks(ticks=[0, 1, 2])
        )
end

end #module
