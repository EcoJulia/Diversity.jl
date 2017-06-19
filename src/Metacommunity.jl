using Compat
using DataFrames

importall Diversity.API

"""
    Metacommunity{FP, A, Part, Sim}

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
type Metacommunity{FP, A, Sim, Part} <: AbstractMetacommunity{FP, A, Sim, Part}
    sourceabundances::A
    internalabundances::A
    types::Sim
    partition::Part
    ordinariness::Nullable{A}

    function
        (::Type{Metacommunity{FP, A,
                              Sim, Part}}){FP <: AbstractFloat,
                                           A <: AbstractArray,
                                           Sim <: AbstractTypes,
                                           Part <: AbstractPartition}(abundances::A,
                                                                      types::Sim,
                                                                      part::Part)
        mcmatch(abundances, types, part) ||
            throw(ErrorException("Type or size mismatch between abundance array, partition and type list"))
        internalabundances = calcabundance(types, abundances)
        new{FP, A, Sim, Part}(abundances, internalabundances,
                              types, part, Nullable{A}())
    end

    function
        (::Type{Metacommunity{FP, A,
                              Sim, Part}}){FP <: AbstractFloat,
                                           A <: AbstractArray,
                                           Sim <: AbstractTypes,
                                           Part <: AbstractPartition}(abundances::A,
                                                                      meta::Metacommunity{FP, A, Sim, Part})
        mcmatch(abundances, meta.types, meta.part) ||
            throw(ErrorException("Type or size mismatch between abundance array, partition and type list"))
        internalabundances = calcabundance(types, abundances)
        new{FP, A, Sim, Part}(abundances, internalabundances, meta.types,
                              meta.part, meta.ordinariness)
    end
end

function Metacommunity{A <: AbstractArray,
                       Meta <: AbstractMetacommunity}(abundances::A, meta::Meta)
    types = gettypes(meta)
    part = getpartition(meta)
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(A), A,
                             typeof(types), typeof(part)}(abundances, types, part)
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
                             typeof(types), typeof(part)}(ab, types, part)
    end
end

function Metacommunity{M <: AbstractMatrix, Sim <: AbstractTypes,
                       Part <: AbstractPartition}(abundances::M,
                                                  types::Sim = UniqueTypes(size(abundances, 1)),
                                                  part::Part = Subcommunities(size(abundances, 2)))
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(M), M,
                             typeof(types), typeof(part)}(abundances, types, part)
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
                             typeof(types), typeof(part)}(ab, types, part)
    end
end

function Metacommunity{V <: AbstractVector, Sim <: AbstractTypes,
                       Part <: AbstractPartition}(abundances::V,
                                                  types::Sim = UniqueTypes(size(abundances, 1)),
                                                  part::Part = Onecommunity())
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(V), V,
                             typeof(types), typeof(part)}(abundances, types, part)
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
                             typeof(types), typeof(part)}(ab, types, part)
    end
end

function Metacommunity{V <: AbstractVector,
                       M <: AbstractMatrix}(abundances::V, zmatrix::M)
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(V), V,
                             GeneralTypes, Onecommunity}(abundances,
                                                         GeneralTypes(zmatrix),
                                                         Onecommunity())
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
                             GeneralTypes,
                             Onecommunity}(ab, GeneralTypes(zmatrix),
                                           Onecommunity())
    end
    
end

function Metacommunity{M <: AbstractMatrix}(abundances::M, zmatrix::M)
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(M), M,
                             GeneralTypes,
                             Subcommunities}(abundances,
                                             GeneralTypes(zmatrix),
                                             Subcommunities(size(abundances, 2)))
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
                             GeneralTypes, Subcommunities}(ab, GeneralTypes(zmatrix),
                                                           Subcommunities(size(ab, 2)))
    end
end

function _gettypes{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.types
end

function _getpartition{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.partition
end

function _getabundance{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part},
                                         input::Bool)
    return input ? meta.sourceabundances : meta.internalabundances
end

function _getordinariness!(meta::Metacommunity)
    if isnull(meta.ordinariness)
        meta.ordinariness = _calcordinariness(meta.types, meta.sourceabundances)
    end
    get(meta.ordinariness)
end
