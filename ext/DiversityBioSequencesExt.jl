# SPDX-License-Identifier: BSD-2-Clause

module DiversityBioSequencesExt

import Diversity
using Diversity.API

using BioSequences
using StringDistances

"""
    GeneticFASTA

Genetic similarity type built from a vector of aligned `BioSequence`s. Each
sequence is a type; similarity is derived from pairwise sequence distances.
"""
struct GeneticFASTA{GeneticData} <: Diversity.AbstractGenetic
    dat::GeneticData
    names::Vector{String}
    ntypes::Int64
    Zmatrix::Matrix{Float64}
end

# Pairwise Hamming distance between aligned sequences.
function _sequencedistances(::Val{:hamming}, dat)
    nseq = length(dat)
    dist = zeros(Float64, nseq, nseq)
    for a in 1:nseq, b in (a + 1):nseq
        dist[a, b] = dist[b, a] = Hamming()(dat[a], dat[b])
    end
    return dist
end

function _sequencedistances(::Val{D}, _) where {D}
    return throw(ArgumentError("unknown sequence distance :$D (try :hamming)"))
end

"""
    GeneticType(dat::AbstractVector{<:BioSequence}; distance = :hamming,
                names = string.(1:length(dat)), transform = :linear,
                k = 1, normalise = true)

Construct a `GeneticFASTA` similarity type from a vector of aligned
`BioSequence`s. Each sequence is a type. `distance` selects the pairwise
sequence distance (`:hamming`) and `transform`, `k` and `normalise` control the
distance-to-similarity conversion (see rdiversity's `dist2sim`).
"""
function Diversity.GeneticType(dat::AbstractVector{S};
                               distance::Symbol = :hamming,
                               names::AbstractVector = string.(1:length(dat)),
                               transform::Symbol = :linear,
                               k::Real = 1,
                               normalise::Bool = true) where {S <: BioSequence}
    dist = _sequencedistances(Val(distance), dat)
    Zmatrix = Diversity._dist2sim(dist; transform = transform, k = k,
                                  normalise = normalise, max_d = maximum(dist))
    return GeneticFASTA{typeof(dat)}(dat, String.(names), length(dat), Zmatrix)
end

import Diversity.API: _getdiversityname
_getdiversityname(::GeneticFASTA) = "Genetic (sequence)"

end
