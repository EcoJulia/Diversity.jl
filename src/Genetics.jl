using Diversity
using Diversity.API

import PopGen
import StringDistances
import LinearAlgebra    # Symmetric()

function _hammingDistance(geno1, geno2)
    ismissing(geno1) || ismissing(geno2) && return missing
    if length(geno1) > 2
        @warn "hamming_distance may not work correctly for ploidy > 2"
    end
    #TODO Fix ploidy > 2 - e.g. (1, 1, 1, 2) ≠ (1, 2, 2, 2)
    
    max(sum(geno1 .∉ Ref(geno2)), sum(geno2 .∉ Ref(geno1)))
end

function genDistance(dat::PopData) 
    # Initialise objects
    matrix_obj = PopGen.loci_matrix(dat)
    N = size(matrix_obj, 1)
    output = zeros(Float64, N, N)
    indices = PopGen.pairwise_pairs(1:N)

    # Calculate distance matrix
    for (a, b) in indices
        output[a, b] = sum(_hammingDistance.((@view matrix_obj[a, :]),
                                             (@view matrix_obj[b, :])))
    end
    return LinearAlgebra.Symmetric(output)
end

function genDistance(dat::AbstractVector) 
    # Initialise objects
    N = length(dat)
    output = zeros(Int64, N, N)
    indices = PopGen.pairwise_pairs(1:N)

    # Calculate distance matrix
    for (a, b) in indices
        output[a, b] = StringDistances.evaluate(StringDistances.Hamming(), dat[a][2], dat[b][2])
    end
    return LinearAlgebra.Symmetric(output)
end