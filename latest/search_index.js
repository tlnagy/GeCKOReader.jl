var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#GeCKOReader.jl-1",
    "page": "Home",
    "title": "GeCKOReader.jl",
    "category": "section",
    "text": "(Image: Project Status)(Image: Build Status)Some analysis code for interpreting GeCKOv2 CRISPR knockout screens"
},

{
    "location": "index.html#GeCKOReader.read_seq_file-Tuple{Dict{Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},String},String,Bio.Seq.BioSequence{Bio.Seq.DNAAlphabet{4}},Int64,Int64,Int64}",
    "page": "Home",
    "title": "GeCKOReader.read_seq_file",
    "category": "Method",
    "text": "read_seq_file(filename, barcode, stagger, read_length; verbose=false)\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.gen_barcode_mapping-Tuple{DataFrames.DataFrame,Int64}",
    "page": "Home",
    "title": "GeCKOReader.gen_barcode_mapping",
    "category": "Method",
    "text": "gen_barcode_mapping(barcodes::DataFrame, len::Int)\n\nGenerate a dictionary mapping of sequence of type Bio.Seq.DNASequence to their string id from barcodes. Truncate all barcodes to length len. The final result is not guaranteed to be unique due to the truncation step.\n\n\n\n"
},

{
    "location": "index.html#Loading-1",
    "page": "Home",
    "title": "Loading",
    "category": "section",
    "text": "These functions are useful for loading fastq.gz filesModules=[GeCKOReader]\nPrivate=false\nPages=[\"io.jl\", \"utils.jl\"]"
},

{
    "location": "index.html#GeCKOReader.calc_log2fc-Tuple{DataFrames.DataFrame,DataFrames.DataFrame}",
    "page": "Home",
    "title": "GeCKOReader.calc_log2fc",
    "category": "Method",
    "text": "calc_log2fc(exp::DataFrame, control::DataFrame)\n\nCombine control with exp and compute the log2 fold change between their :rel_freqs columns\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.compute_gene_stats-Tuple{DataFrames.DataFrame,Array{Float64,1}}",
    "page": "Home",
    "title": "GeCKOReader.compute_gene_stats",
    "category": "Method",
    "text": "compute_gene_stats(df::DataFrame, negcontrols_fc::Vector{Float64};\n                   column::Symbol=:log2fc)\n\nCalculates pvalues and medians on a per gene basis for the given DataFrame df using the specified column. P-values are computed using a Mann-Whitney test of the gene values against the values provided in negcontrols_fc for the negative controls.\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.fit_lines-Tuple{DataFrames.DataFrame,DataFrames.DataFrame,DataFrames.DataFrame}",
    "page": "Home",
    "title": "GeCKOReader.fit_lines",
    "category": "Method",
    "text": "fit_lines(t0::DataFrame, t1::DataFrame, t2::DataFrame)\n\nFit a line to the relative frequencies of each guide at t0, t1, and t2 by solving the linear least squares equation for an overdetermined system.\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.get_log2fc_pvals-Tuple{DataFrames.DataFrame,DataFrames.DataFrame}",
    "page": "Home",
    "title": "GeCKOReader.get_log2fc_pvals",
    "category": "Method",
    "text": "get_log2fc_pvals(exp::DataFrame, control::DataFrame)\n\nGiven a treatment dataframe exp and a control dataframe control calculates the log 2 fold changes on a per guide level and then uses a Mann Whitney test to compare this the negative control guides to determine the significance on a per gene level.\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.get_rel_freqs-Tuple{GeCKOReader.ScreenSample,DataFrames.DataFrame}",
    "page": "Home",
    "title": "GeCKOReader.get_rel_freqs",
    "category": "Method",
    "text": "get_rel_freqs(sample::ScreenSample,\n                   callable_bc::DataFrame)::DataFrame\n\nComputes the relative frequencies of barcodes as compared to the median frequency of occurence for all the negative controls. The frequencies are computed initially by dividing by the total number of reads for that sample.\n\n\n\n"
},

{
    "location": "index.html#Processing-1",
    "page": "Home",
    "title": "Processing",
    "category": "section",
    "text": "Modules=[GeCKOReader]\nPrivate=false\nPages=[\"processing.jl\"]"
},

{
    "location": "index.html#GeCKOReader.Plot.plot_all_qualities-Tuple{GeCKOReader.ScreenSample}",
    "page": "Home",
    "title": "GeCKOReader.Plot.plot_all_qualities",
    "category": "Method",
    "text": "plot_all_qualities(s::ScreenSample)\n\nGenerate a 2D histogram of read qualities from the FASTQ file per position in the read using all of the read info.\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.Plot.plot_guide_trajectories-Tuple{DataFrames.DataFrame,String}",
    "page": "Home",
    "title": "GeCKOReader.Plot.plot_guide_trajectories",
    "category": "Method",
    "text": "plot_guide_trajectories(lines::DataFrame, gene::String)\n\nPlot the trajectories of the guides for a specific gene given the data in lines. The lines dataframe should be formatted like the output of fit_lines.\n\n\n\n"
},

{
    "location": "index.html#GeCKOReader.Plot.plot_hit_qualities-Tuple{GeCKOReader.ScreenSample}",
    "page": "Home",
    "title": "GeCKOReader.Plot.plot_hit_qualities",
    "category": "Method",
    "text": "plot_hit_qualities(s::ScreenSample)\n\nGenerate a 2D histogram of read qualities from the FASTQ file per position in the read using only reads that were matched.\n\n\n\n"
},

{
    "location": "index.html#Plotting-1",
    "page": "Home",
    "title": "Plotting",
    "category": "section",
    "text": "Modules=[GeCKOReader.Plot]\nPages=[\"plotting.jl\"]"
},

]}
