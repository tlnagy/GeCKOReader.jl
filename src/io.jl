"""
read_seq_file(filename, barcode, stagger, read_length; verbose=false)

"""
function read_seq_file(mapping::Dict{DNASequence, String},
                       filename::String,
                       barcode::DNASequence,
                       stagger::Int,
                       read_length::Int,
                       truncate::Int;
                       verbose=false)
    # inflate streamed data from zipped file
    stream = open(filename) |> ZlibInflateInputStream
    try
        reader = FASTQReader(stream, quality_encoding=:illumina18)
        record = FASTQSeqRecord{DNASequence}()

        searchseq = barcode * dna"TCTTGTGGAAAGGACGAAACACCG"
        len_search = length(searchseq)
        # for detecting alternative staggers
        search_query = ExactSearchQuery(searchseq)

        sample = ScreenSample(filename, read_length, 41)

        while !eof(reader)
            read!(reader, record)
            sample.reads += 1
            # number of errors in the constant region
            dist = hamming(searchseq, record.seq[(stagger+1):(stagger+len_search)])
            # get quality scores
            qual = record.metadata.quality

            priming_range = search(record.seq, search_query)
            (priming_range.stop > 0) && (sample.alt_staggers[priming_range.start] += 1)

            if verbose && sample.reads < 100
                println("\n\n", record.seq, " Read")
                println("*"^stagger, searchseq, "*"^(read_length - len_search - stagger), " Target constant region")
                println("-"^stagger, record.seq[stagger+1:stagger+len_search], "-"^(read_length - len_search - stagger), " Read constant region")
                println("-"^(stagger+len_search), record.seq[stagger+len_search+1:end], " Read sgrna")
                println(length(record.seq[stagger+len_search+1:end]))
            end

            if dist <= 5
                # read must have a PHRED quality above 37 for entire barcode region
                sgrna_quality = qual[(stagger + len_search+1):(stagger+len_search + truncate)]

                if sum(sgrna_quality .< 37) == 0
                    for pos in 1:read_length
                        sample.hit_qualities[qual[pos], pos] += 1
                    end

                    sgrna = record.seq[(stagger + len_search + 1):(stagger + len_search + truncate)]

                    if haskey(mapping, sgrna)
                        sample.hits += 1
                        (verbose && sample.reads < 100) && println("\nMatched as $(mapping[sgrna])")
                        if haskey(sample.counts, mapping[sgrna])
                            sample.counts[mapping[sgrna]] += 1
                        else
                            sample.counts[mapping[sgrna]] = 1
                        end

                    else
                        (verbose && sample.reads < 100) && println("\nNo match for $sgrna")
                        sample.no_match += 1
                    end
                elseif verbose && sample.reads < 100
                    println("\nDoes not pass quality threshold in sgrna\n$(join(sgrna_quality, ','))")
                else
                    sample.low_quality += 1
                end
            elseif verbose && sample.reads < 100
                println("\nToo many errors in constant region ($(dist) errors detected)")
            else
                sample.constant_error += 1
            end

            for pos in 1:read_length
                sample.all_qualities[qual[pos], pos] += 1
            end
        end
        # compare the provided stagger versus the most common detected one and insure that they
        # are the same. this is a sanity check that we're not throwing away too many reads.
        alt_stagger = findmax(sample.alt_staggers)[2] - 1
        (alt_stagger != stagger) && warn("An alternative stagger of $alt_stagger was detected as the dominant stagger versus the provided stagger of $stagger.")
        return sample
    finally
        close(stream)
        end
end
