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
    compute_gene_stats(df::DataFrame, negcontrols_fc::Vector{Float64};
                       column::Symbol=:log2fc)

Calculates pvalues and medians on a per gene basis for the given DataFrame `df`
using the specified `column`. P-values are computed using a Mann-Whitney test of
the gene values against the values provided in `negcontrols_fc` for the negative
controls.
"""
function compute_gene_stats(df::DataFrame, negcontrols_fc::Vector{Float64};
                            column::Symbol=:log2fc)
    res = @from i in df begin
        @where isfinite(_try_get(getfield(i, column)))
        @group i by i.gene into g
        @let test = -log10(pvalue(MannWhitneyUTest(map(j->_try_get(getfield(j, column)), g),
                                                   negcontrols_fc)))
        @let med = median(map(j->_try_get(getfield(j, column)), g))
        @select {
            gene = _try_get(g.key),
            med_val = med,
            pval = test,
            prod = abs(med * test),
            n_guides = length(g)
        }
        @collect DataFrame
    end
    rename!(res, :med_val, Symbol("med_", column))
    res
end


"""
    fit_lines(t0::DataFrame, t1::DataFrame, t2::DataFrame)

Fit a line to the relative frequencies of each guide at `t0`, `t1`, and `t2`
by solving the linear least squares equation for an overdetermined system.
"""
function fit_lines(t0::DataFrame, t1::DataFrame, t2::DataFrame)
    data = join(join(t0, t1, on=:id), t2, on=:id)

    X = hcat(ones(3, 1), [0,1,2])
    tmp = @from i in data begin
        # Solve the system of linear equations for the betahats, which give you the slope
        # and the intercept
        # https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)#The_general_problem
        @let betahat = X \ [get(i.rel_freqs), get(i.rel_freqs_1), get(i.rel_freqs_2)]
        @select {
            id=get(i.id),
            gene=get(i.gene),
            intercept=betahat[1],
            slope=betahat[2],
            t0_count=get(i.counts),
            t2_count=get(i.counts_2)
        }
        @collect DataFrame
    end
end
