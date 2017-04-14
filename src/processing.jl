"""
    get_rel_freqs(sample::ScreenSample,
                       callable_bc::DataFrame)::DataFrame

Computes the relative frequencies of barcodes as compared to the median
frequency of occurence for all the negative controls. The frequencies are
computed initially by dividing by the total number of reads for that sample.
"""
function get_rel_freqs(sample::ScreenSample,
                       callable_bc::DataFrame)::DataFrame
    # construct a DataFrame from the dictionary
    df = DataFrame(id=collect(keys(sample.counts)),
                   counts=float(collect(values(sample.counts))),
                   sample_id=sample.filename)

    total_reads = sum(df[:counts])


    # use all rows from barcodes even if they are missing
    tmp = join(callable_bc, df, on=:id, kind=:left)
    # set missing sgRNAs to the pseudocount value
    tmp[isna(tmp[:counts]), :counts] = 0
    # compute frequencies by normalizing by the total number of reads
    tmp[:freqs] = tmp[:counts] ./ total_reads

    negcontrols = @from i in tmp begin
        @where contains(get(i.gene), "NonTargetingControlGuide") && i.freqs > 0
        @select i
        @collect DataFrame
    end
    (size(negcontrols)[1] < 100) && warn("Very few negative controls detected, this might be problem")

    negcontrol_med = median(negcontrols[:freqs])

    tmp[:rel_freqs] = tmp[:freqs] ./ negcontrol_med

    tmp
end

"""
    get_log2fc_pvals(exp::DataFrame, control::DataFrame)

Given a treatment dataframe `exp` and a control dataframe `control` calculates
the log 2 fold changes on a per guide level and then uses a Mann Whitney test
to compare this the negative control guides to determine the significance on
a per gene level.
"""
function get_log2fc_pvals(exp::DataFrame, control::DataFrame)

    # join controls dataframe with treatment and compute log2fc between them
    x = calc_log2fc(exp, control)

    # get all log2fc values for the negative controls
    negcontrols_fc = @from i in x begin
        @where contains(get(i.gene), "NonTargetingControlGuide") && isfinite(get(i.log2fc))
        @select get(i.log2fc)
        @collect
    end

    z = compute_gene_stats(x, negcontrols_fc)

    sort!(z, cols=:pval, rev=true)
    z[:ranking] = collect(1:size(z)[1])
    z

end

"""
    calc_log2fc(exp::DataFrame, control::DataFrame)

Combine `control` with `exp` and compute the log2 fold change between their
`:rel_freqs` columns
"""
function calc_log2fc(exp::DataFrame, control::DataFrame)
    @from i in control begin
        @join j in exp on i.id equals j.id
        @select {
                    i.id,
                    log2fc=Nullable(log2(get(j.rel_freqs)/get(i.rel_freqs))),
                    rel_freq1=i.rel_freqs,
                    rel_freq2=j.rel_freqs,
                    i.gene,
                    i.sample_id
                }
        @collect DataFrame
    end
end


"""
    compute_gene_stats(df::DataFrame, negcontrols_fc::Vector{Float64})

Calculates pvalues and median log2 fold changes on a per gene basis for the
given DataFrame `df` and the vector of log2 fold changes for the negative
controls.
"""
function compute_gene_stats(df::DataFrame, negcontrols_fc::Vector{Float64})
    @from i in df begin
        @where isfinite(get(i.log2fc))
        @group i by i.gene into g
        @let test = -log10(pvalue(MannWhitneyUTest(map(j->get(j.log2fc), g), negcontrols_fc)))
        @let med = median(map(j->get(j.log2fc), g))
        @select {
            gene = get(g.key),
            med_log2fc = med,
            pval = test,
            prod = abs(med * test),
            n_guides = length(g)
        }
        @collect DataFrame
    end
end
