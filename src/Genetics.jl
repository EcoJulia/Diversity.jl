using Diversity
using Diversity.API

using PopGen
using StringDistances
using LinearAlgebra

using FASTX

abstract type AbstractGeneticTypes{PopData} <: 
    Diversity.API.AbstractTypes
end

struct GeneticType{PopData} <: AbstractGeneticTypes{PopData}
    dat::PopData
    ntypes::Int64
    Zmatrix::Matrix{Float64}
end

function _hammingDistance(geno1, geno2)
    ismissing(geno1) || ismissing(geno2) && return missing
    if length(geno1) > 2
        @warn "hamming_distance may not work correctly for ploidy > 2"
    end
    #TODO Fix ploidy > 2 - e.g. (1, 1, 1, 2) ≠ (1, 2, 2, 2)
    
    max(sum(geno1 .∉ Ref(geno2)), sum(geno2 .∉ Ref(geno1)))
end

function GeneticType(dat::PopData) 
    # Initialise objects
    matrix_obj = PopGen.loci_matrix(dat)
    ntypes = size(matrix_obj, 1)
    output = zeros(Float64, ntypes, ntypes)
    indices = PopGen.pairwise_pairs(1:ntypes)

    # Calculate distance matrix
    for (a, b) in indices
        output[a, b] = sum(_hammingDistance.((@view matrix_obj[a, :]),
                                             (@view matrix_obj[b, :])))
    end
    dist = Symmetric(output)
    dist /= maximum(dist)

    # Calculate similarity matrix
    Zmatrix = 1 .- dist

    return GeneticType{PopData}(dat, ntypes, Zmatrix)
end

function GeneticType(dat::Vector) # Vector{BioSequences.AminoAcidSequence}
    # Initialise objects
    ntypes = length(dat)
    output = zeros(Int64, ntypes, ntypes)
    indices = PopGen.pairwise_pairs(1:ntypes)

    # Calculate distance matrix
    for (a, b) in indices
        output[a, b] = evaluate(Hamming(), dat[a], dat[b])
    end
    dist = Symmetric(output)
    dist /= maximum(dist)

    # Calculate similarity matrix
    Zmatrix = 1 .- dist

    return GeneticType{BioSequences.AminoAcidSequence}(dat, ntypes, Zmatrix)
end