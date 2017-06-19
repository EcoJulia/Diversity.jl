using Compat
using Diversity

### AbstractPartition API

"""
    AbstractPartition

Abstract supertype for all partitioning types. AbstractPartition
subtypes allow you to define how to partition your total metacommunity
(e.g. an ecosystem) into smaller components (e.g. subcommunities).
"""
@compat abstract type AbstractPartition end

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
function _countsubcommunities(p::AbstractPartition)
    return length(_getsubcommunitynames(p))
end

### AbstractTypes API

"""
    AbstractTypes

Abstract supertype for all similarity types. Its subtypes allow you to
define how similarity is measured between individuals.
"""
@compat abstract type AbstractTypes end

"""
    _gettypenames(t::AbstractTypes, input::Bool)

Returns the names of the types in an AbstractTypes object. Must be
implemented by each AbstractTypes subtype. `input` determines whether
to count the number of input or output types, which varies, for
instance, when the types are determined by a phylogeny.
"""
function _gettypenames end

"""
    _calcsimilarity(t::AbstractTypes, a::AbstractArray)

Retrieves (and possibly calculates) a similarity matrix from t. Must be
implemented by each AbstractTypes subtype.
"""
function _calcsimilarity end

"""
    _counttypes(::AbstractTypes, input::Bool)

Returns number of types in an AbstractTypes object, t. May be
implemented by each AbstractTypes subtype. `input` determines whether
to count the number of input or output types, which varies, for
instance, when the types are determined by a phylogeny. Default is to
count length of corresponding types name vector.
"""
function _counttypes end

function _counttypes(t::AbstractTypes, input::Bool)
    return length(_gettypenames(t, input))
end

"""
    _calcabundance(t::AbstractTypes, a::AbstractArray)

Calculates the abundance a for AbstractTypes, t (if necessary). May be
implemented by each AbstractTypes subtype.
"""
function _calcabundance end
function _calcabundance(::AbstractTypes, a::AbstractArray)
    return a
end

"""
    _calcordinariness(t::AbstractTypes, a::AbstractArray)

Calculates the ordinariness of abundance a from AbstractTypes, t. May be
implemented by each AbstractTypes subtype.
"""
function _calcordinariness end
function _calcordinariness(t::AbstractTypes, a::AbstractArray)
    return _calcsimilarity(t, a) * _calcabundance(t, a)
end

### AbstractMetacommunity API

"""
    AbstractMetacommunity{FP, A, Sim, Part}

AbstractMetacommunity is the abstract supertype of all metacommunity
types. AbstractMetacommunity subtypes allow you to define how to
partition your total metacommunity (e.g. an ecosystem) into smaller
components (e.g. subcommunities), and how to assess similarity between
individuals within it.

"""
@compat abstract type AbstractMetacommunity{FP <: AbstractFloat,
    A <: AbstractArray,
    Sim <: AbstractTypes,
    Part <: AbstractPartition} end

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
    _getabundance(m::AbstractMetacommunity, input::Bool)

Returns the abundances array of the metacommunity. Must be implemented
by each AbstractMetacommunity subtype.
"""
function _getabundance end

"""
    _getmetaabundance(m::AbstractMetacommunity, input::Bool)

Returns the metacommunity abundances of the metacommunity. May be
implemented by each AbstractMetacommunity subtype.
"""
function _getmetaabundance end
function _getmetaabundance(m::AbstractMetacommunity, input::Bool)
    return reduce(+, SubcommunityIterator(m, mc -> _getabundance(mc, input)))
end

"""
    _getweight(m::AbstractMetacommunity)

Returns the subcommunity weights of the metacommunity. May be
implemented by each AbstractMetacommunity subtype.
"""
function _getweight end
function _getweight(m::AbstractMetacommunity)
    return reduce(+, TypeIterator(m, mc -> _getabundance(mc, false)))'
end

"""
    _getordinariness!(m::AbstractMetacommunity, a)

Returns (and possibly calculates) the ordinariness array of the
subcommunities. May be implemented by each AbstractMetacommunity
subtype.
"""
function _getordinariness! end
function _getordinariness!(m::AbstractMetacommunity)
    return _calcordinariness(_gettypes(m), _getabundance(m, false))
end

"""
    _getmetaordinariness!(m::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness of the
metacommunity as a whole. May be implemented by each
AbstractMetacommunity subtype.
"""
function _getmetaordinariness! end
function _getmetaordinariness!(m::AbstractMetacommunity)
    return reduce(+, SubcommunityIterator(m, _getordinariness!))
end

### Other optional APIs to implement

"""
    floattypes(t)

This function returns a set containing the floating point types that
are compatible with the Diversity-related object, t.

"""
function floattypes end

function floattypes{A <: AbstractArray}(::A)
    return Set([eltype(A)])
end

function floattypes(::AbstractTypes)
    return Set(subtypes(AbstractFloat))
end

function floattypes(::AbstractPartition)
    return Set(subtypes(AbstractFloat))
end

function floattypes{FP, A, Sim, Part}(::AbstractMetacommunity{FP, A, Sim, Part})
    return Set([FP])
end

"""
    typematch(args...)

Checks whether the types of a variety of Diversity-related objects
have compatible types (using floattypes()).

"""
typematch(args...) = length(mapreduce(floattypes, ∩, args)) ≥ 1

"""
    mcmatch(a::AbstractArray, part::AbstractPartition, sim::AbstractTypes)

Checks for type and size compatibility for elements contributing to a Metacommunity
"""
function mcmatch end

function mcmatch(m::AbstractMatrix, sim::AbstractTypes, part::AbstractPartition)
    realm = _calcabundance(sim, m)
    return typematch(realm, sim, part) &&
        _counttypes(sim, false) == size(realm, 1) &&
        _countsubcommunities(part) == size(realm, 2) &&
        sum(realm) ≈ 1
end

function mcmatch(v::AbstractVector, sim::AbstractTypes, part::AbstractPartition)
    realv = _calcabundance(sim, v)
    return typematch(realv, sim, part) &&
        _counttypes(sim, false) == size(realv, 1) &&
        _countsubcommunities(part) == 1 &&
        sum(realv) ≈ 1
end

"""
    vectorise(arg)

Returns the argument if it is an vector, or reduce an array to a
vector, or return a vector containing the argument if it's a number

"""
function vectorise end

@inline function vectorise(arr::AbstractVector)
  arr
end

@inline function vectorise(arr::AbstractArray)
  vec(arr)
end

@inline function vectorise(num::Real)
  [num]
end
