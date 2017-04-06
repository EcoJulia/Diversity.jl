using Compat

"""
    floattypes(t)

This function returns a set containing the floating point types that
are compatible with the Diversity-related object, t.

"""
function floattypes end

function floattypes{A <: AbstractArray}(::A)
    return Set([eltype(A)])
end

"""
    typematch(args...)

Checks whether the types of a variety of Diversity-related objects
have compatible types (using floattypes()).

"""
typematch(args...) = length(mapreduce(floattypes, ∩, args)) ≥ 1

"""
    AbstractPartition

Abstract supertype for all partitioning types. AbstractPartition
subtypes allow you to define how to partition your total metacommunity
(e.g. an ecosystem) into smaller components (e.g. subcommunities).

"""
@compat abstract type AbstractPartition end

"""
    countsubcommunities(p::AbstractPartition)

Returns number of subcommunities in a partition, p.

"""
function countsubcommunities end

"""
    getnames(arg)

Returns the names of the subcommunities of the AbstractPartition or
the names of the types of the AbstractTypes.

"""
function getnames end

function getnames(::AbstractPartition)
    throw(NullException())
end

function floattypes(::AbstractPartition)
    return Set(subtypes(AbstractFloat))
end

"""
    AbstractTypes

Abstract supertype for all similarity types. Its subtypes allow you to
define how similarity is measured between individuals.

"""
@compat abstract type AbstractTypes end

"""
    counttypes(t::AbstractTypes)

Returns number of types in an AbstractTypes object, t.

"""
function counttypes end

"""
    getsimilarity(t::AbstractTypes)

Retrieves (and possibly calculates) a similarity matrix from p

"""
function getsimilarity end

"""
    getordinariness(t::AbstractTypes, a::AbstractArray)

Calculates the ordinariness of abundance a from AbstractTypes, t

"""
function getordinariness end

function getordinariness(t::AbstractTypes, a::AbstractArray)
    getsimilarity(t) * a
end

function getnames(::AbstractTypes)
    throw(NullException())
end

function floattypes(::AbstractTypes)
    return Set(subtypes(AbstractFloat))
end

"""
    mcmatch(ab::AbstractArray, part::AbstractPartition, sim::AbstractTypes)

Checks for type and size compatibility for elements contributing to a Metacommunity
"""
function mcmatch end

function mcmatch(ab::AbstractMatrix, sim::AbstractTypes, part::AbstractPartition)
    typematch(ab, sim, part) &&
    counttypes(sim) == size(ab, 1) &&
    countsubcommunities(part) == size(ab, 2) &&
    sum(ab) ≈ 1
end

function mcmatch(ab::AbstractVector, sim::AbstractTypes, part::AbstractPartition)
    typematch(ab, sim, part) &&
    counttypes(sim) == size(ab, 1) &&
    countsubcommunities(part) == 1 &&
    sum(ab) ≈ 1
end

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
    getabundance(::AbstractMetacommunity)

Returns the abundances array of the metacommunity.

"""
function getabundance end

"""
    getordinariness!(::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness array of the subcommunities.

"""
function getordinariness! end

"""
    getmetaordinariness!(::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness of the metacommunity as a whole.

"""
function getmetaordinariness! end

@inline function getmetaordinariness!{Meta <: AbstractMetacommunity}(meta::Meta)
    sumoversubcommunities(meta, getordinariness!(meta))
end

"""
    getweight(::AbstractMetacommunity)

Retrieves (and possibly calculates) the relative weights of the subcommunities.

"""
@inline function getweight{Meta <: AbstractMetacommunity}(meta::Meta)
    sumovertypes(meta, getabundance(meta))
end

"""
    getpartition(::AbstractMetacommunity)

Returns the AbstractPartition component of the metacommunity.

"""
function getpartition end

"""
    gettypes(::AbstractMetacommunity)

Returns the AbstractTypes component of the metacommunity.

"""
function gettypes end

function floattypes{FP, A, Sim, Part}(::AbstractMetacommunity{FP, A, Sim, Part})
    return Set([FP])
end

"""
    sumovertypes(::AbstractMetacommunity)

Sums an array over its types, leaving an array of the same
dimensionality (but length of 1st dimension is 1).

"""
@inline function sumovertypes{FP, A, Sim, Part}(meta::AbstractMetacommunity{FP, A, Sim, Part}, values::A)
    mapslices(sum, values, 1)::A
end

"""
    sumoversubcommunities(::AbstractMetacommunity)

Sums an array over its subcommunities, leaving an array of the same
dimensionality (but length of 2nd and subsequent dimensions are 1).

"""
function sumoversubcommunities end

@inline function sumoversubcommunities{FP, A, Sim, Part}(meta::AbstractMetacommunity{FP, A, Sim, Part}, values::A)
    mapslices(sum, values, collect(2:ndims(values)))::A
end

@inline function sumoversubcommunities{FP, V <: AbstractVector, Sim, Part}(meta::AbstractMetacommunity{FP, V, Sim, Part}, values::V)
    values::V
end


## Now create the functions for the iterator interface for Partitions
## and Metacommunities
import Base.start, Base.next, Base.done, Base.eltype, Base.length

@inline function start(meta::AbstractMetacommunity)
    1::Int64
end

@inline function next{Meta <: AbstractMetacommunity}(meta::Meta, state::Int64)
    (view(getabundance(meta), :, state), state + 1)
end

@inline function done{Meta <: AbstractMetacommunity}(meta::Meta, state::Int64)
    index > countsubcommunities(getpartition(meta))
end

@inline function eltype{FP, A, Sim, Part}(::AbstractMetacommunity{FP, A, Sim, Part})
    FP
end

@inline function length{Meta <: AbstractMetacommunity}(meta::Meta)
    countsubcommunities(getpartition(meta))
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
