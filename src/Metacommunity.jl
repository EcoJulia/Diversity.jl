# SPDX-License-Identifier: BSD-2-Clause

using DataFrames
using Missings
using EcoBase

"""
    Metacommunity{FP, ARaw, AProcessed, Sim, Part}

Metacommunity type, representing a whole metacommunity containing a
single community or a collection of subcommunities. The metacommunity
of individuals *may* be further partitioned into smaller groups. For
instance this may be an ecosystem, which consists of a series of
subcommunities. The AbstractPartition subtype within it stores
relative abundances of different types, e.g. species, and also allows
for similarity between individuals.

# Constructor:

Metacommunity(abundances::AbstractArray,
              types::AbstractTypes,
              part::AbstractPartition)

# Members:

- `abundances` the abundance matrix for the metacommunity.

- `partition` the instance of the AbstractPartition subtype, containing the
  subcommunities.

- `types` The instance of the AbstractTypes subtype, from which
  similarities between individuals can be calculated.

- `ordinariness` A cache of the ordinariness of the individuals in the
  Partition. Should only be accessed through
  getordinariness!(::Metacommunity), which will populate the cache if
  it has not yet been calculated.

- `weights` A cache of the subcommunity weights, accessed through
  getweight(::Metacommunity).

- `metaordinariness` A cache of the ordinariness of the metacommunity
  as a whole, accessed through getmetaordinariness!(::Metacommunity).

All three caches are populated on first use and never invalidated, so a
metacommunity should not be mutated once it has been measured.
"""
mutable struct Metacommunity{FP, ARaw, AProcessed, Sim, Part} <:
               Diversity.API.AbstractMetacommunity{FP, ARaw, AProcessed, Sim,
                                                   Part}
    rawabundances::ARaw
    processedabundances::AProcessed
    scale::FP
    types::Sim
    partition::Part
    ordinariness::Union{AProcessed, Missing}
    weights::Union{Vector{FP}, Missing}
    metaordinariness::Union{Vector{FP}, Missing}

    function Metacommunity{FP, ARaw, AProcessed,
                           Sim, Part}(abundances::ARaw,
                                      matrix::AProcessed,
                                      types::Sim,
                                      part::Part) where
        {FP <: AbstractFloat, ARaw <: AbstractArray,
         AProcessed <: AbstractArray{FP},
         Sim <: AbstractTypes, Part <: AbstractPartition}
        mcmatch(matrix, types, part) ||
            error("Type or size mismatch between abundance array, " *
                  "partition and type list")
        processedabundances, scale = _calcabundance(types, matrix)
        return new{FP, ARaw, AProcessed, Sim, Part}(abundances,
                                                    processedabundances, scale,
                                                    types, part, missing,
                                                    missing, missing)
    end
end

function Metacommunity(abundances::ARaw,
                       meta::Meta) where
    {ARaw <: AbstractArray, Meta <: AbstractMetacommunity}
    types = gettypes(meta)
    part = getpartition(meta)
    mat = reshape(abundances, counttypes(types), countsubcommunities(part))
    if sum(mat) ≉ one(eltype(mat))
        @warn "Abundances not normalised to 1, correcting..."
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), ARaw, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types,
                                                      part)
end

function Metacommunity(abundances::V,
                       types::Sim = UniqueTypes(size(abundances, 1)),
                       part::Part = Onecommunity()) where
    {V <: AbstractVector, Sim <: AbstractTypes, Part <: AbstractPartition}
    mat = reshape(abundances / sum(abundances), length(abundances), 1)
    return Metacommunity{eltype(mat), V, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types,
                                                      part)
end

function Metacommunity(abundances::V,
                       types::Sim = UniqueTypes(size(abundances, 1)),
                       part::Part = Onecommunity()) where
    {FP <: AbstractFloat, V <: AbstractVector{FP},
     Sim <: AbstractTypes, Part <: AbstractPartition}
    mat = reshape(abundances, length(abundances), 1)
    if sum(mat) ≉ one(eltype(mat))
        @warn "Abundances not normalised to 1, correcting..."
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), typeof(abundances), typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types,
                                                      part)
end

function Metacommunity(abundances::M,
                       types::Sim = UniqueTypes(size(abundances, 1)),
                       part::Part = Subcommunities(size(abundances, 2))) where
    {M <: AbstractMatrix, Sim <: AbstractTypes, Part <: AbstractPartition}
    mat = abundances / sum(abundances)
    return Metacommunity{eltype(mat), M, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types,
                                                      part)
end

function Metacommunity(abundances::M,
                       types::Sim = UniqueTypes(size(abundances, 1)),
                       part::Part = Subcommunities(size(abundances, 2))) where
    {FP <: AbstractFloat, M <: AbstractMatrix{FP},
     Sim <: AbstractTypes, Part <: AbstractPartition}
    mat = abundances
    if sum(mat) ≉ one(eltype(mat))
        @warn "Abundances not normalised to 1, correcting..."
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), typeof(abundances), typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types,
                                                      part)
end

function Metacommunity(abundances::V,
                       zmatrix::M) where
    {FP <: AbstractFloat, V <: AbstractVector, M <: AbstractMatrix{FP}}
    return Metacommunity(abundances, GeneralTypes(zmatrix), Onecommunity())
end

function Metacommunity(abundances::MU,
                       zmatrix::M) where
    {FP <: AbstractFloat, MU <: AbstractMatrix, M <: AbstractMatrix{FP}}
    return Metacommunity(abundances, GeneralTypes(zmatrix),
                         Subcommunities(size(abundances, 2)))
end

# Keep a partition that is already one of ours. A foreign `AbstractPlaces` is not an
# `AbstractPartition` and cannot be, so rebuild one from its names instead — `placenames` is part of
# EcoBase's own interface for `AbstractPlaces`, so it is always there to ask, and `collect` is what
# turns whatever vector of strings it returns into the `Vector{String}` `Subcommunities` takes.
_aspartition(part::AbstractPartition) = part
function _aspartition(places::EcoBase.AbstractPlaces)
    return Subcommunities(collect(String, placenames(places)))
end

function Metacommunity(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) || return Metacommunity(occurrences(asm))

    # Materialise the similarity as a plain matrix and hand it to GeneralTypes, whatever the
    # original types were, so the result computes similarity-sensitive diversity through the
    # ordinary path with no dependence on the source hierarchy. The *processed* abundances go with
    # it: for a phylogeny those are the branch abundances the scaled Zmatrix is indexed by, and
    # `calcsimilarity(t, getscale(asm))` is exactly the matrix the measures would have used.
    processed = getabundance(asm)
    types = GeneralTypes(calcsimilarity(gettypes(asm), _getscale(asm)),
                         gettypenames(asm))
    return Metacommunity(processed, types, _aspartition(getpartition(asm)))
end

import Diversity.API._gettypes
_gettypes(meta::Metacommunity) = meta.types

import Diversity.API._getpartition
_getpartition(meta::Metacommunity) = meta.partition

import Diversity.API._getabundance
function _getabundance(meta::Metacommunity, raw::Bool)
    return raw ? meta.rawabundances : meta.processedabundances
end

import Diversity.API._getordinariness!
function _getordinariness!(meta::Metacommunity)
    # Bound to a local before the test and returned from the local, rather than re-read from the
    # field: the field is declared `Union{AProcessed, Missing}`, so returning it directly infers as
    # that union, while a local is narrowed by `ismissing` to the array alone.
    ord = meta.ordinariness
    if ismissing(ord)
        ord = _calcordinariness(meta.types, meta.processedabundances,
                                meta.scale)
        meta.ordinariness = ord
    end
    return ord
end

# The subcommunity weights and the metacommunity ordinariness are cached for the same reason the
# ordinariness itself is: each is a full reduction over an ntypes x nplaces array -- 4.3 ms and
# 4.9 ms respectively at 200 types x 200,000 places -- and *every* measure built over the
# metacommunity asks for them, so `diversity` over the seven measures repeated both seven times.
# Only the matrix case is cached here; where the raw abundances are a vector there is one
# subcommunity, so the weight is `[1]` and the metacommunity ordinariness is the subcommunity
# ordinariness, already cached. Those fall through to the defaults in `API.jl`.
import Diversity.API._getweight
function _getweight(meta::Metacommunity{FP, <:AbstractMatrix}) where {FP}
    w = meta.weights
    if ismissing(w)
        summed = sum(_getabundance(meta, false), dims = 1)
        w = reshape(summed, length(summed))
        meta.weights = w
    end
    return w
end

import Diversity.API._getmetaordinariness!
function _getmetaordinariness!(meta::Metacommunity{FP,
                                                   <:AbstractMatrix}) where {FP}
    mord = meta.metaordinariness
    if ismissing(mord)
        summed = sum(_getordinariness!(meta), dims = 2)
        mord = reshape(summed, length(summed))
        meta.metaordinariness = mord
    end
    return mord
end

import Diversity.API._getscale
_getscale(m::Metacommunity) = m.scale
