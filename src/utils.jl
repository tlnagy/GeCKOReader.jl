"""
    hamming(a::DNASequence, b::DNASequence)

Compute number of mismatches between two DNA sequences of equal length
"""
function Distances.hamming(a::DNASequence, b::DNASequence)
    n = length(a)
    (n != length(b)) && throw(DimensionMismatch("Lengths of sequences must be equal"))
    mismatches = 0
    for i in 1:n
        if a[i] != b[i]
            mismatches += 1
        end
    end
    mismatches
end


"""
    gen_barcode_mapping(barcodes::DataFrame, len::Int)

Generate a dictionary mapping of sequence of type Bio.Seq.DNASequence to their
string id from `barcodes`. Truncate all barcodes to length `len`. The final
result is not guaranteed to be unique due to the truncation step.
"""
function gen_barcode_mapping(barcodes::DataFrame, len::Int)
    mapping = Dict{DNASequence, String}()
    for row in eachrow(barcodes)
        mapping[DNASequence(row[:sequence][1:len])] = row[:id]
    end
    mapping
end

_try_get(x) = x
_try_get(x::Nullable) = get(x)
