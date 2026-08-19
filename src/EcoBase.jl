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
import EcoBase: view
using EcoBase: asindices
occurrences(mc::AbstractMetacommunity) = getabundance(mc)
places(mc::AbstractMetacommunity) = getpartition(mc)
things(mc::AbstractMetacommunity) = gettypes(mc)

# What the units are called, so that anything printing a metacommunity says what it means rather than
# "thing" and "place". EcoBase supplies these four hooks and defaults them on the *assemblage*; we
# answer on the **types** and the **partition** instead, because the unit is a property of what is
# being measured, not of the metacommunity that holds it. That is what lets `PhyloBranches` say
# "branch" — its things really are branches, and the extension overrides these there.
import EcoBase: thingkind, thingkindplural, placekind, placekindplural

thingkind(mc::AbstractMetacommunity) = thingkind(gettypes(mc))
thingkind(::AbstractTypes) = "species"
thingkindplural(mc::AbstractMetacommunity) = thingkindplural(gettypes(mc))
thingkindplural(::AbstractTypes) = "species"

placekind(mc::AbstractMetacommunity) = placekind(getpartition(mc))
placekind(::AbstractPartition) = "subcommunity"
placekindplural(mc::AbstractMetacommunity) = placekindplural(getpartition(mc))
placekindplural(::AbstractPartition) = "subcommunities"

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

# Subsetting, which is what `view` below is made of.
#
# The default materialises the similarity matrix and hands it to a GeneralTypes. That is right for
# any types object whose `_calcsimilarity` is a matrix, which is all of them except `UniqueTypes`,
# whose identity similarity is a `UniformScaling` and cannot be indexed -- it supplies its own
# method. Note the scale is baked into the matrix, so the result needs no scale of its own: that is
# exactly what makes a subset of a phylogeny's branches go on measuring correctly, and it is why an
# arbitrary subset of branches becomes a GeneralTypes rather than staying a phylogeny. It is not a
# tree any more.
import Diversity.API: _subsettypes, _subsetpartition

function _subsettypes(t::AbstractTypes, idx, scale::Real)
    return GeneralTypes(calcsimilarity(t, scale)[idx, idx],
                        gettypenames(t, false)[idx])
end

function _subsetpartition(p::AbstractPartition, idx)
    return Subcommunities(getsubcommunitynames(p)[idx])
end

"""
    SubAssemblage(types, partition, occurrences)

A subset of a metacommunity, as returned by `view`. It is an EcoBase
`AbstractAssemblage` rather than a `Metacommunity`, which is what lets it be a
genuine view: it aliases the parent's abundances rather than copying them, and
it carries them unnormalised, since only a `Metacommunity` requires them to sum
to one. The diversity measures still work on it, because abundances are
normalised when they are read.
"""
struct SubAssemblage{FP <: Real, T <: AbstractTypes, P <: AbstractPartition,
                     A <: AbstractArray{FP}} <:
       EcoBase.AbstractAssemblage{FP, T, P}
    types::T
    partition::P
    occurrences::A
    scale::Float64
end

things(sub::SubAssemblage) = sub.types
_getscale(sub::SubAssemblage) = sub.scale
places(sub::SubAssemblage) = sub.partition
occurrences(sub::SubAssemblage) = sub.occurrences

# The units are the subset's types', not "thing" and "place" -- see the naming hooks above.
Base.show(io::IO, sub::SubAssemblage) = _showassemblage(io, sub)

thingkind(sub::SubAssemblage) = thingkind(things(sub))
thingkindplural(sub::SubAssemblage) = thingkindplural(things(sub))
placekind(sub::SubAssemblage) = placekind(places(sub))
placekindplural(sub::SubAssemblage) = placekindplural(places(sub))

# The scale a subset should measure with. It is its own, not the parent's: for a phylogeny the scale
# is the abundance-weighted mean root-to-tip distance, so dropping subcommunities changes it. It
# cannot be recovered from the branch abundances alone, which is why it is worked out here, from the
# leaf abundances of the subcommunities kept, while the parent metacommunity is still to hand.
# Everything except a phylogeny has a scale of one, and short-circuits.
function _subsetscale(mc::AbstractMetacommunity, st)
    _getscale(mc) == 1 && return 1.0
    raw = getabundance(mc, true)
    raw isa AbstractMatrix || return Float64(_getscale(mc))
    kept = raw[:, st]
    return Float64(_calcabundance(gettypes(mc), kept ./ sum(kept))[2])
end

# Whether a selection keeps everything, in order -- in which case there is nothing to subset.
function _keepsall(idx, n)
    return length(idx) == n && all(i == j for (i, j) in zip(idx, Base.OneTo(n)))
end

"""
    view(mc::AbstractMetacommunity; species, sites)

Takes a view of part of a metacommunity, returning a `SubAssemblage` over the
types at `species` and the subcommunities at `sites`. Either may be given as
indices, as a boolean mask, or as names, and either may be omitted to keep
everything.

Warning: `species` selects **things**, whatever the metacommunity's types call
them. For a phylogeny a thing is a *branch*, not a species, because
`PhyloBranches` measures over branches -- ask `EcoBase.thingkind` if in doubt, or
look at what `show` calls them. The keyword name is EcoBase's.

A dimension you do not restrict keeps its object untouched, so restricting only
`sites` leaves a phylogeny a phylogeny. Restricting `species` cannot: an
arbitrary subset of branches is no longer a tree, so the types become a
`GeneralTypes` carrying the scaled similarity submatrix, which measures
identically but reports itself as `"Arbitrary Z"` and drops any added output
columns.

The result aliases the parent's abundances rather than copying them, so it
reflects later changes to them, and unlike a `Metacommunity` it caches nothing.
Its abundances are a subset of the parent's and so do not sum to one; the
measures normalise them on reading, which means the subset measures as a
metacommunity in its own right. Use `Metacommunity(view(...))` for a converted,
cached object instead.
"""
function view(mc::AbstractMetacommunity;
              species = Base.OneTo(counttypes(mc)),
              sites = Base.OneTo(countsubcommunities(mc)))
    # `asindices` is EcoBase's own name resolution, so names, masks and indices all work here
    # exactly as they do for the rest of its interface.
    sp = asindices(species, gettypenames(mc))
    st = asindices(sites, getsubcommunitynames(mc))
    # A dimension that is not actually restricted keeps its object untouched, which is not just an
    # optimisation: subsetting types is what turns a phylogeny into a GeneralTypes, so a view that
    # only picks subcommunities must not do it, or restricting sites would silently cost you the
    # tree.
    scale = _subsetscale(mc, st)
    # When the types are subset the scale is baked into the similarity matrix, so the result needs
    # none of its own; when they are passed through it has to be carried instead.
    keeptypes = _keepsall(sp, counttypes(mc))
    types = keeptypes ? gettypes(mc) : _subsettypes(gettypes(mc), sp, scale)
    part = _keepsall(st, countsubcommunities(mc)) ? getpartition(mc) :
           _subsetpartition(getpartition(mc), st)
    return SubAssemblage(types, part, Base.view(occurrences(mc), sp, st),
                         keeptypes ? scale : 1.0)
end

RecipesBase.@recipe function f(var::DataFrame, asm::AbstractAssemblage)
    return var[!, :diversity], getcoords(places(asm))
end
