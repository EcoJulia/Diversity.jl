using Compat
using DataFrames

importall Diversity.API

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
type Metacommunity{FP, ARaw, AProcessed, Sim, Part} <:
    AbstractMetacommunity{FP, ARaw, AProcessed, Sim, Part}
    rawabundances::ARaw
    processedabundances::AProcessed
    scale::FP
    types::Sim
    partition::Part
    ordinariness::Nullable{AProcessed}

    function
        (::Type{Metacommunity{FP, ARaw, AProcessed, Sim,
                              Part}}){FP <: AbstractFloat,
                                      ARaw <: AbstractArray,
                                      AProcessed <: AbstractArray,
                                      Sim <: AbstractTypes,
                                      Part <: AbstractPartition}(abundances::ARaw,
                                                                 matrix::AProcessed,
                                                                 types::Sim,
                                                                 part::Part)
        mcmatch(matrix, types, part) ||
            error("Type or size mismatch between abundance array, partition and type list")
        processedabundances, scale = _calcabundance(types, matrix)
        new{FP, ARaw, AProcessed, Sim, Part}(abundances, processedabundances, scale,
                                             types, part, Nullable{AProcessed}())
    end
end

function Metacommunity{A <: AbstractArray,
                       Meta <: AbstractMetacommunity}(abundances::A, meta::Meta)
    types = gettypes(meta)
    part = getpartition(meta)
    mat = reshape(abundances, counttypes(types), countsubcommunities(part))
    if sum(mat) ≉ one(eltype(mat))
        warn("Abundances not normalised to 1, correcting...")
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), A, typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types, part)
end

function Metacommunity{Sim <: AbstractTypes,
                       Part <: AbstractPartition}(abundances::Vector,
                                                  types::Sim =
                                                  UniqueTypes(size(abundances, 1)),
                                                  part::Part =
                                                  Onecommunity())
    mat = reshape(abundances, length(abundances), 1)
    if sum(mat) ≉ one(eltype(mat))
        warn("Abundances not normalised to 1, correcting...")
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), typeof(abundances), typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types, part)
end

function Metacommunity{M <: AbstractMatrix, Sim <: AbstractTypes,
                       Part <: AbstractPartition}(abundances::M,
                                                  types::Sim =
                                                  UniqueTypes(size(abundances, 1)),
                                                  part::Part =
                                                  Subcommunities(size(abundances, 2)))
    mat = abundances
    if sum(mat) ≉ one(eltype(mat))
        warn("Abundances not normalised to 1, correcting...")
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), typeof(abundances), typeof(mat),
                         typeof(types), typeof(part)}(abundances, mat, types, part)
end

function Metacommunity{M <: AbstractMatrix}(abundances::Vector, zmatrix::M)
    mat = reshape(abundances, length(abundances), 1)
    if sum(mat) ≉ one(eltype(mat))
        warn("Abundances not normalised to 1, correcting...")
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), typeof(abundances), typeof(mat),
                         GeneralTypes, Onecommunity}(abundances, mat,
                                                     GeneralTypes(zmatrix),
                                                     Onecommunity())
end

function Metacommunity{M <: AbstractMatrix}(abundances::M, zmatrix::M)
    mat = abundances
    if sum(mat) ≉ one(eltype(mat))
        warn("Abundances not normalised to 1, correcting...")
        mat = mat / sum(mat)
    end
    return Metacommunity{eltype(mat), typeof(abundances), typeof(mat),
                         GeneralTypes,
                         Subcommunities}(abundances, mat,
                                         GeneralTypes(zmatrix),
                                         Subcommunities(size(abundances, 2)))
end

function _gettypes{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.types
end

function _getpartition{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.partition
end

function _getabundance{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part},
                                         raw::Bool)
    return raw ? meta.rawabundances : meta.processedabundances
end

function _getordinariness!(meta::Metacommunity)
    if isnull(meta.ordinariness)
        meta.ordinariness = _calcordinariness(meta.types, meta.processedabundances, meta.scale)
    end
    get(meta.ordinariness)
end

function _getscale(m::Metacommunity)
    return m.scale
end
