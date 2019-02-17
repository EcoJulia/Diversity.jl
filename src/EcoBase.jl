using Diversity.API
using EcoBase
using EcoBase: AbstractAssemblage, AbstractThings, AbstractPlaces
using Compat
using RecipesBase
using DataFrames

# Now satisfy the EcoBase interface
import EcoBase: nplaces, placenames
nplaces(part::AbstractPartition) = countsubcommunities(part)
placenames(part::AbstractPartition) = getsubcommunitynames(part)

import EcoBase: coordinates
# Where there's no spatial information, make it up!
function coordinates(part::AbstractPartition{Nothing})
    dimx = round(Int, sqrt(countsubcommunities(part)), RoundUp)
    coords = Matrix{Float64}(undef, countsubcommunities(part), 2)
    for i in Base.OneTo(countsubcommunities(part))
        coords[i, 1] = (i - 1) % dimx + 1
        coords[i, 2] = (i - 1) ÷ dimx + 1
    end
    return coords
end

import EcoBase: nthings, thingnames
nthings(types::AbstractTypes) = counttypes(types)
thingnames(types::AbstractTypes) = gettypenames(types)

import EcoBase: occurrences, places, things
occurrences(mc::AbstractMetacommunity) = getabundance(mc)
places(mc::AbstractMetacommunity) = getpartition(mc)
things(mc::AbstractMetacommunity) = gettypes(mc)

# And use the EcoBase interface to provide a basic diversity interface
import Diversity.API: _getpartition
_getpartition(p::AbstractAssemblage) = places(p)

import Diversity.API: _gettypes
_gettypes(p::AbstractAssemblage) = things(p)

import Diversity.API: _getaddedoutput
_getaddedoutput(::AbstractThings) = nothing

import Diversity.API: _addedoutputcols
_addedoutputcols(::AbstractThings) = Dict{Symbol, Type}()

import Diversity.API: _getscale
_getscale(::AbstractAssemblage) = 1

import Diversity.API: _getabundance
_getabundance(a::AbstractAssemblage, raw::Bool) =
    raw ? occurrences(a) : occurrences(a) / sum(occurrences(a))

import Diversity.API: _getsubcommunitynames
_getsubcommunitynames(p::AbstractPlaces) = placenames(p)

import Diversity.API: _countsubcommunities
_countsubcommunities(p::AbstractPlaces) = length(_getsubcommunitynames(p))

import Diversity.API: _gettypenames
_gettypenames(p::AbstractThings, ::Bool) = thingnames(p)

import Diversity.API: _counttypes
_counttypes(t::AbstractThings, raw::Bool) = length(_gettypenames(t, raw))

import Diversity.API: _calcsimilarity
_calcsimilarity(t::AbstractThings, ::Real) =
    Matrix(1.0I, counttypes(t), counttypes(t))

import Diversity.API: _getweight
function _getweight(a::AbstractAssemblage)
    ab = _getabundance(a, false)
    w = Compat.sum(ab, dims=1)
    return reshape(w, length(w))
end

import Diversity.API: _calcabundance
function _calcabundance(::AbstractThings, a::AbstractArray)
    return a, one(eltype(a))
end

import Diversity.API: _calcordinariness
function _calcordinariness(t::AbstractThings, a::AbstractArray, ::Real)
    abundance, scale = _calcabundance(t, a)
    return _calcsimilarity(t, scale) * abundance
end

import Diversity.API: _getordinariness!
function _getordinariness!(a::AbstractAssemblage)
    return _calcordinariness(_gettypes(a), _getabundance(a, false),
                             _getscale(a))
end

import Diversity.API: _getmetaordinariness!
function _getmetaordinariness!(a::AbstractAssemblage)
    ord = Compat.sum(_getordinariness!(a), dims=2)
    return reshape(ord, length(ord))
end

import Diversity.API: _getdiversityname
_getdiversityname(::AbstractThings) = "species"

import Diversity.API: _hassimilarity
_hassimilarity(::AbstractThings) = false

RecipesBase.@recipe function f(var::DataFrame, asm::AbstractAssemblage)
    var[:diversity], getcoords(places(asm))
end
