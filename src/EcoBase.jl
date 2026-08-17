# SPDX-License-Identifier: BSD-2-Clause

using Diversity.API
using EcoBase
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

import EcoBase: getcoords
# EcoBase's fallback is `getcoords(plc::AbstractPlaces{Nothing}) = plc`, on the grounds that a
# partition with no location data has no coordinates to give. Ours *does* — `coordinates` above makes
# up a grid — so without this method the plot recipes hand Plots a partition object instead of a
# coordinate matrix, and plotting any metacommunity without real spatial data fails outright.
getcoords(part::AbstractPartition{Nothing}) = coordinates(part)

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
function _getabundance(a::AbstractAssemblage, raw::Bool)
    return raw ? occurrences(a) : occurrences(a) / sum(occurrences(a))
end

import Diversity.API: _getsubcommunitynames
_getsubcommunitynames(p::AbstractPlaces) = placenames(p)

import Diversity.API: _countsubcommunities
_countsubcommunities(p::AbstractPlaces) = length(_getsubcommunitynames(p))

import Diversity.API: _gettypenames
_gettypenames(p::AbstractThings, ::Bool) = thingnames(p)

import Diversity.API: _counttypes
_counttypes(t::AbstractThings, raw::Bool) = length(_gettypenames(t, raw))

import Diversity.API: _calcsimilarity
# The similarity of types that have none: everything is like itself and nothing else. Shared with
# the AbstractTypes guard at the end of the file, which hands back the same matrix.
_identitysimilarity(t) = Matrix(1.0I, counttypes(t), counttypes(t))

function _calcsimilarity(t::AbstractThings, ::Real)
    return _identitysimilarity(t)
end

import Diversity.API: _getweight
function _getweight(a::AbstractAssemblage)
    ab = _getabundance(a, false)
    w = sum(ab, dims = 1)
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
    ord = sum(_getordinariness!(a), dims = 2)
    return reshape(ord, length(ord))
end

import Diversity.API: _getdiversityname
_getdiversityname(::AbstractThings) = "species"

import Diversity.API: _hassimilarity
_hassimilarity(::AbstractThings) = false

# The reverse bridge above is typed on EcoBase's supertypes, and Diversity's own abstract types are
# subtypes of them — so without these guards a Diversity subtype that implements *neither* side
# recurses (the forward method calls the underscore API, which dispatches straight back to the
# reverse method) until the stack overflows. These more specific methods break the cycle and say
# what is actually missing instead. They must stay: any subtype that *does* implement the underscore
# API defines a still more specific method and wins, so these are only ever reached by an incomplete
# implementation.

# Reports the API function an incomplete implementation failed to provide, in place of whatever the
# EcoBase fallback would otherwise have done silently — `why` says what that was.
function _notimplemented(fname, x,
                         why = "the EcoBase fallback cannot be used here " *
                               "because it would call back into this method")
    return error("$fname is not implemented for $(typeof(x)). A Diversity " *
                 "subtype must implement it; $why.")
end

_getpartition(m::AbstractMetacommunity) = _notimplemented("_getpartition", m)
_gettypes(m::AbstractMetacommunity) = _notimplemented("_gettypes", m)
function _getabundance(m::AbstractMetacommunity, ::Bool)
    return _notimplemented("_getabundance", m)
end
function _getsubcommunitynames(p::AbstractPartition)
    return _notimplemented("_getsubcommunitynames", p)
end
_gettypenames(t::AbstractTypes, ::Bool) = _notimplemented("_gettypenames", t)

# `_calcsimilarity` fails differently, and needs its own guard: the AbstractThings fallback above
# does not recurse, it quietly returns an identity matrix. That is *right* for a type declaring it
# has no similarity — which is why the fallback exists — but for one that claims similarity and
# never supplied the matrix it means silently wrong answers instead of an error. So make the claim,
# not the omission, the thing that is refused. `_hassimilarity` is a trait with one method per type,
# so the branch below is folded away during inference rather than tested at run time.
function _calcsimilarity(t::AbstractTypes, ::Real)
    _hassimilarity(t) &&
        _notimplemented("_calcsimilarity", t,
                        "its `_hassimilarity` reports that it has similarity, " *
                        "so the identity matrix the EcoBase fallback would " *
                        "return is not right for it")
    return _identitysimilarity(t)
end

RecipesBase.@recipe function f(var::DataFrame, asm::AbstractAssemblage)
    return var[!, :diversity], getcoords(places(asm))
end
