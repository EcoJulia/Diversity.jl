using DataFrames
using Missings
using Compat: @warn
using EcoBase

"""
    Metacommunity{FP, ARaw, AProcessed, Part, Sim}

Metacommunity type, representing a whole metacommunity containing a
single community or a collection of subcommunities. The metacommunity
of individuals *may* be further partitioned into smaller groups. For
instance this may be an ecosystem, which consists of a series of
subcommunities. The AbstractPartition subtype within it stores
relative abundances of different types, e.g. species, and also allows
for similarity between individuals.

# Constructor:

Metacommunity(abundances::AbstractArray,
              part::AbstractPartition,
              types::AbstractTypes)

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

"""
mutable struct Metacommunity{FP, ARaw, AProcessed, Sim, Part} <:
    Diversity.API.AbstractMetacommunity{FP, ARaw, AProcessed, Sim, Part}
    rawabundances::ARaw
    processedabundances::AProcessed
    scale::FP
    types::Sim
    partition::Part
    ordinariness::Union{AProcessed, Missing}

    function Metacommunity{FP, ARaw, AProcessed,
                           Sim, Part}(abundances::ARaw,
                                      matrix::AProcessed,
                                      types::Sim,
                                      part::Part) where
        {FP <: AbstractFloat, ARaw <: AbstractArray, AProcessed <: AbstractArray{FP},
         Sim <: AbstractTypes, Part <: AbstractPartition}
        mcmatch(matrix, types, part) ||
            error("Type or size mismatch between abundance array, " *
                  "partition and type list")
        processedabundances, scale = _calcabundance(types, matrix)
        new{FP, ARaw, AProcessed, Sim, Part}(abundances, processedabundances, scale,
                                             types, part, missing)
    end
end

function Metacommunity(abundances::ARaw, meta::Meta) where
    {ARaw <: AbstractArray, Meta <: AbstractMetacommunity}
    types = gettypes(meta)
    part = getpartition(meta)
    mat = reshape(abundances, counttypes(types), countsubcommunities(part))
    if sum(mat) ≉ one(eltype(mat))
        @warn "Abundances not normalised to 1, correcting..."
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), ARaw, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types, part)
end

function Metacommunity(abundances::V,
                       types::Sim = UniqueTypes(size(abundances, 1)),
                       part::Part = Onecommunity()) where
    {V <: AbstractVector, Sim <: AbstractTypes, Part <: AbstractPartition}
    mat = reshape(abundances / sum(abundances), length(abundances), 1)
    return Metacommunity{eltype(mat), V, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types, part)
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
                         typeof(types), typeof(part)}(abundances, mat, types, part)
end

function Metacommunity(abundances::M,
                       types::Sim = UniqueTypes(size(abundances, 1)),
                       part::Part = Subcommunities(size(abundances, 2))) where
    {M <: AbstractMatrix, Sim <: AbstractTypes, Part <: AbstractPartition}
    mat = abundances / sum(abundances)
    return Metacommunity{eltype(mat), M, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types, part)
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
                         typeof(types), typeof(part)}(abundances, mat, types, part)
end

function Metacommunity(abundances::V, zmatrix::M) where
    {FP <: AbstractFloat, V <: AbstractVector, M <: AbstractMatrix{FP}}
    return Metacommunity(abundances, GeneralTypes(zmatrix), Onecommunity())
end

function Metacommunity(abundances::MU, zmatrix::M) where
    {FP <: AbstractFloat, MU <: AbstractMatrix{FP}, M <: AbstractMatrix{FP}}
    return Metacommunity(abundances, GeneralTypes(zmatrix),
                         Subcommunities(size(abundances, 2)))
end

Metacommunity(asm::EcoBase.AbstractAssemblage) =
    hassimilarity(asm) ?
        Metacommunity(occurrences(asm), _calcsimilarity(asm)) :
        Metacommunity(occurrences(asm))

import Diversity.API._gettypes
function _gettypes(meta::Metacommunity{FP, ARaw, AProcessed, Sim, Part}) where
    {FP, ARaw, AProcessed, Sim, Part}
    return meta.types
end

import Diversity.API._getpartition
function _getpartition(meta::Metacommunity{FP, ARaw, AProcessed, Sim, Part}) where
    {FP, ARaw, AProcessed, Sim, Part}
    return meta.partition
end

import Diversity.API._getabundance
function _getabundance(meta::Metacommunity{FP, ARaw, AProcessed, Sim, Part},
                       raw::Bool) where {FP, ARaw, AProcessed, Sim, Part}
    return raw ? meta.rawabundances : meta.processedabundances
end

import Diversity.API._getordinariness!
function _getordinariness!(meta::Metacommunity)
    if ismissing(meta.ordinariness)
        meta.ordinariness = _calcordinariness(meta.types, meta.processedabundances, meta.scale)
    end
    return meta.ordinariness
end

import Diversity.API._getscale
function _getscale(m::Metacommunity)
    return m.scale
end
