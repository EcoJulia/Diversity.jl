using Diversity.API

"""
    gettypes(m::AbstractMetacommunity)

Returns the AbstractTypes component of the metacommunity.
"""
function gettypes(m::AbstractMetacommunity)
    return _gettypes(m)
end

"""
    getpartition(m::AbstractMetacommunity)

Returns the AbstractPartition component of the metacommunity.

"""
function getpartition(m::AbstractMetacommunity)
        return _getpartition(m)
end

"""
    counttypes(m::AbstractMetacommunity[, raw::Bool = false])
    counttypes(t::AbstractTypes[, raw::Bool = false])

Returns number of types in an `AbstractTypes` object or the
`AbstractMetacommunity` containing it. `raw` determines whether to
count the number of raw or processed types, which varies, for instance,
when the types are determined by a phylogeny.
"""
function counttypes end

function counttypes(m::AbstractMetacommunity, raw::Bool = false)
    return _counttypes(_gettypes(m), raw)
end

function counttypes(t::AbstractTypes, raw::Bool = false)
    return _counttypes(t, raw)
end

"""
    countsubcommunities(m::AbstractMetacommunity)
    countsubcommunities(p::AbstractPartition)

Returns number of subcommunities in an `AbstractPartition` object or the
`AbstractMetacommunity` containing it.
"""
function countsubcommunities end

function countsubcommunities(m::AbstractMetacommunity)
    return _countsubcommunities(_getpartition(m))
end

function countsubcommunities(p::AbstractPartition)
    return _countsubcommunities(p)
end

"""
    gettypenames(m::AbstractMetacommunity[, raw::Bool = false])
    gettypenames(t::AbstractTypes[, raw::Bool = false])

Returns the names of the types of the `AbstractTypes` object or the
`AbstractMetacommunity` containing it. `raw` determines whether to
count the number of raw or processed types, which varies, for instance,
when the types are determined by a phylogeny.
"""
function gettypenames end

function gettypenames(m::AbstractMetacommunity, raw::Bool = false)
    return _gettypenames(_gettypes(m), raw)
end

function gettypenames(t::AbstractTypes, raw::Bool = false)
    return _gettypenames(t, raw)
end

"""
    getdiversityname(m::AbstractMetacommunity)
    getdiversityname(t::AbstractTypes)

Returns the name of the diversity type used.
"""
function getdiversityname end

getdiversityname(m::AbstractMetacommunity) = _getdiversityname(_gettypes(m))
getdiversityname(t::AbstractTypes) = _getdiversityname(t)

"""
    getsubcommunitynames(m::AbstractMetacommunity)
    getsubcommunitynames(p::AbstractPartition)

Returns the names of the subcommunities in an `AbstractPartition` object or the
`AbstractMetacommunity` containing it.
"""
function getsubcommunitynames end

function getsubcommunitynames(m::AbstractMetacommunity)
    return _getsubcommunitynames(_getpartition(m))
end

function getsubcommunitynames(p::AbstractPartition)
    return _getsubcommunitynames(p)
end

"""
    getabundance(m::AbstractMetacommunity)

Returns the abundances array of the metacommunity.
"""
function getabundance(m::AbstractMetacommunity, raw::Bool = false)
    return _getabundance(m, raw)
end

"""
    getmetaabundance(m::AbstractMetacommunity)

Returns the metacommunity abundances of the metacommunity.
"""
function getmetaabundance(m::AbstractMetacommunity, raw::Bool = false)
    return _getmetaabundance(m, raw)
end

"""
    getweight(m::AbstractMetacommunity)

Returns the subcommunity weights of the metacommunity.
"""
function getweight(m::AbstractMetacommunity)
    return _getweight(m)
end

"""
    getordinariness!(m::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness array of the subcommunities.
"""
function getordinariness!(m::AbstractMetacommunity)
    return _getordinariness!(m)
end

"""
    getmetaordinariness!(m::AbstractMetacommunity)

Returns (and possibly calculates) the ordinariness of the
metacommunity as a whole.
"""
function getmetaordinariness!(m::AbstractMetacommunity)
    return _getmetaordinariness!(m)
end

"""
    calcsimilarity(t::AbstractTypes, scale::Real)

Retrieves (and possibly calculates) a similarity matrix from t.
"""
function calcsimilarity(t::AbstractTypes, scale::Real)
    return _calcsimilarity(t, scale)
end

# Now satisfy the EcoBase interface
import Diversity.EcoBase: occurrences, places, things
occurrences(mc::AbstractMetacommunity) = getabundance(mc)
places(mc::AbstractMetacommunity) = getpartition(mc)
things(mc::AbstractMetacommunity) = gettypes(mc)

import Diversity.EcoBase: nplaces, placenames
nplaces(part::AbstractPartition) = countsubcommunities(part)
placenames(part::AbstractPartition) = getsubcommunitynames(part)

import Diversity.EcoBase: nthings, thingnames
nthings(types::AbstractTypes) = counttypes(types)
thingnames(types::AbstractTypes) = gettypenames(types)
