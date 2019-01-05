using Diversity.API
using EcoBase

# Now satisfy the EcoBase interface
import EcoBase: nplaces, placenames, getcoords
nplaces(part::AbstractPartition) = countsubcommunities(part)
placenames(part::AbstractPartition) = getsubcommunitynames(part)
function getcoords(part::AbstractPartition)
    coords = Matrix{Float64}(undef, countsubcommunities(part), 2)
    coords[:, 1] = Base.OneTo(countsubcommunities(part))
    coords[:, 2] .= 1.0
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
import Diversity.API: _getscale
_getscale(::EcoBase.AbstractAssemblage) = 1

import Diversity.API: _getabundance
_getabundance(a::EcoBase.AbstractAssemblage, raw::Bool) =
    raw ? occurrences(a) : occurrences(a) / sum(occurrences)

import Diversity.API: _getsubcommunitynames
_getsubcommunitynames(p::EcoBase.AbstractPlaces) = placenames(p)

import Diversity.API: _gettypenames
_gettypenames(p::EcoBase.AbstractThings) = typenames(p)

import Diversity.API: _calcsimilarity
_calcsimilarity(t::EcoBase.AbstractThings, ::Real) =
    Matrix(1.0I, _counttypes(t), _counttypes(t))
