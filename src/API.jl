using Diversity
using Compat.InteractiveUtils
using Compat
using EcoBase

### Abstract supertype definitions
"""
    AbstractPartition

Abstract supertype for all partitioning types. AbstractPartition
subtypes allow you to define how to partition your total metacommunity
(e.g. an ecosystem) into smaller components (e.g. subcommunities).
"""
abstract type AbstractPartition{LT <: Union{Nothing,
                                            EcoBase.AbstractLocationData}} <:
    EcoBase.AbstractPlaces{LT}
end

"""
    AbstractTypes

Abstract supertype for all similarity types. Its subtypes allow you to
define how similarity is measured between individuals.
"""
abstract type AbstractTypes <: EcoBase.AbstractThings end

"""
    AbstractMetacommunity{FP <: AbstractFloat,
                          ARaw <: AbstractArray,
                          AProcessed <: AbstractMatrix{FP},
                          Sim <: AbstractTypes,
                          Part <: AbstractPartition}

AbstractMetacommunity is the abstract supertype of all metacommunity
types. AbstractMetacommunity subtypes allow you to define how to
partition your total metacommunity (e.g. an ecosystem) into smaller
components (e.g. subcommunities), and how to assess similarity between
individuals within it.

"""
abstract type AbstractMetacommunity{FP <: AbstractFloat,
                                    ARaw <: AbstractArray,
                                    AProcessed <: AbstractMatrix{FP},
                                    Sim <: AbstractTypes,
                                    Part <: AbstractPartition} <:
    EcoBase.AbstractAssemblage{FP, Sim, Part}
end

### AbstractPartition API
"""
    _getsubcommunitynames(p::AbstractPartition)

Returns the names of the subcommunities in the partition
object. Must be implemented by each AbstractPartition subtype.
"""
function _getsubcommunitynames end

"""
    _countsubcommunities(::AbstractPartition)

Returns number of subcommunities in a partition, p. May be implemented
by each AbstractPartition subtype. Default is to count length of
subcommunity name vector.
"""
function _countsubcommunities end

### AbstractTypes API
"""
    _gettypenames(t::AbstractTypes, raw::Bool)

Returns the names of the types in an AbstractTypes object. Must be
implemented by each AbstractTypes subtype. `raw` determines whether
to count the number of raw or processed types, which varies, for
instance, when the types are determined by a phylogeny.
"""
function _gettypenames end

"""
    _calcsimilarity(t::AbstractTypes, scale::Real)

Retrieves (and possibly calculates) a similarity matrix from t. Must be
implemented by each AbstractTypes subtype.
"""
function _calcsimilarity end

"""
    _counttypes(::AbstractTypes, raw::Bool)

Returns number of types in an AbstractTypes object, t. May be
implemented by each AbstractTypes subtype. `raw` determines whether
to count the number of raw or processed types, which varies, for
instance, when the types are determined by a phylogeny. Default is to
count length of corresponding types name vector.
"""
function _counttypes end

"""
    _getdiversityname(::AbstractTypes)

Returns the name of the diversity type used to calculate.
"""
function _getdiversityname end
_getdiversityname(::AbstractTypes) = "unknown"

"""
    _addedoutputcols(::AbstractTypes)

Returns the name of any additional columns needed to be added to outputs.
"""
function _addedoutputcols end

"""
    _getaddedoutput(::AbstractTypes)

Returns the name of any additional columns needed to be added to outputs.
"""
function _getaddedoutput end

"""
    _calcabundance(t::AbstractTypes, a::AbstractArray)

Calculates the abundance a for AbstractTypes, t (if necessary). May be
implemented by each AbstractTypes subtype.
"""
function _calcabundance end
function _calcabundance(::T, a::A) where {T <: AbstractTypes,
                                          A <: AbstractArray}
    return a, one(eltype(a))
end

"""
    _calcordinariness(t::AbstractTypes, a::AbstractArray, scale::Real)

Calculates the ordinariness of abundance a from AbstractTypes, t. May be
implemented by each AbstractTypes subtype.
"""
function _calcordinariness end
function _calcordinariness(t::T, a::A, ::Real) where {T <: AbstractTypes,
                                                      A <: AbstractArray}
    abundance, scale = _calcabundance(t, a)
    return _calcsimilarity(t, scale) * abundance
end

### AbstractMetacommunity API
"""
    _gettypes(::AbstractMetacommunity)

Returns the AbstractTypes component of the metacommunity. Must be
implemented by each AbstractMetacommunity subtype.
"""
function _gettypes end

"""
    _getpartition(::AbstractMetacommunity)

Returns the AbstractPartition component of the metacommunity. Must be
implemented by each AbstractMetacommunity subtype.
"""
function _getpartition end

"""
    _getabundance(m::AbstractMetacommunity, raw::Bool)

Returns the abundances array of the metacommunity. Must be implemented
by each AbstractMetacommunity subtype.
"""
function _getabundance end

"""
    _getmetaabundance(m::AbstractMetacommunity, raw::Bool)

Returns the metacommunity abundances of the metacommunity. May be
implemented by each AbstractMetacommunity subtype.
"""
function _getmetaabundance end
_getmetaabundance(mc::Meta, raw::Bool) where
{FP, AProcessed, Sim, Part,
 Meta <: Diversity.API.AbstractMetacommunity{FP, <: AbstractVector,
                                             AProcessed, Sim, Part}} =
    _getabundance(mc, raw)

function _getmetaabundance(mc::Meta, raw::Bool) where
    {FP, AProcessed, Sim, Part,
     Meta <: Diversity.API.AbstractMetacommunity{FP, <: AbstractMatrix,
                                                 AProcessed, Sim, Part}}
    ab = Compat.sum(_getabundance(mc, raw), dims=2)
    return reshape(ab, length(ab))
end

"""
    _getscale(m::AbstractMetacommunity)

Returns a scaling factor for the metacommunity (needed for
phylogenetics). Normally ignored. Must be implemented by each
AbstractMetacommunity subtype.
"""
function _getscale end

"""
    _getweight(m::AbstractMetacommunity)

Returns the subcommunity weights of the metacommunity. May be
implemented by each AbstractMetacommunity subtype.
"""
function _getweight end
_getweight(::Meta) where
{FP, AProcessed, Sim, Part,
 Meta <: Diversity.API.AbstractMetacommunity{FP, <: AbstractVector,
                                             AProcessed, Sim, Part}} = [one(FP)]

function _getweight(mc::Meta) where
    {FP, AProcessed, Sim, Part,
     Meta <: Diversity.API.AbstractMetacommunity{FP, <: AbstractMatrix,
                                                 AProcessed, Sim, Part}}
    ab = _getabundance(mc, false)
    w = Compat.sum(ab, dims=1)
    return reshape(w, length(w))
end

"""
    _getordinariness!(m::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness array of the
subcommunities. May be implemented by each AbstractMetacommunity
subtype.
"""
function _getordinariness! end

"""
    _getmetaordinariness!(m::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness of the
metacommunity as a whole. May be implemented by each
AbstractMetacommunity subtype.
"""
function _getmetaordinariness! end
_getmetaordinariness!(mc::Meta) where
{FP, AProcessed, Sim, Part,
 Meta <: Diversity.API.AbstractMetacommunity{FP, <: AbstractVector,
                                             AProcessed, Sim, Part}} =
    _getordinariness!(mc)

function _getmetaordinariness!(mc::Meta) where
    {FP, AProcessed, Sim, Part,
     Meta <: Diversity.API.AbstractMetacommunity{FP, <: AbstractMatrix,
                                                 AProcessed, Sim, Part}}
    ord = Compat.sum(_getordinariness!(mc), dims=2)
    return reshape(ord, length(ord))
end

### Other optional APIs to implement

_hassimilarity(::Diversity.API.AbstractTypes) = true

"""
    floattypes(t)

This function returns a set containing the floating point types that
are compatible with the Diversity-related object, t.

"""
function floattypes end

function floattypes(::A) where {FP <: AbstractFloat, A <: AbstractArray{FP}}
    return Set([FP])
end

function floattypes(::T) where T <: AbstractTypes
    return Set(subtypes(AbstractFloat))
end

function floattypes(::P) where P <: AbstractPartition
    return Set(subtypes(AbstractFloat))
end

function floattypes(::M) where
    {FP, ARaw, AProcessed, Sim, Part,
     M <: AbstractMetacommunity{FP, ARaw, AProcessed, Sim, Part}}
    return Set([FP])
end

"""
    typematch(args...)

Checks whether the types of a variety of Diversity-related objects
have compatible types (using floattypes()).

"""
typematch(args...) = length(mapreduce(floattypes, ∩, args)) ≥ 1

"""
    mcmatch(procm::AbstractArray, sim::AbstractTypes, part::AbstractPartition)

Checks for type and size compatibility for elements contributing to a Metacommunity
"""
function mcmatch end

function mcmatch(procm::M, sim::T, part::P) where {M <: AbstractMatrix,
                                                   T <: AbstractTypes,
                                                   P <: AbstractPartition}
    realm = _calcabundance(sim, procm)[1]
    return typematch(realm, sim, part) &&
        _counttypes(sim, true) == size(procm, 1) &&
        _counttypes(sim, false) == size(realm, 1) &&
        _countsubcommunities(part) == size(realm, 2) &&
        sum(realm) ≈ 1.0
end
