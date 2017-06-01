using Compat
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
    counttypes(t::AbstractTypes)

Returns number of types in an AbstractTypes object, t.
"""
function counttypes(t::AbstractTypes, input::Bool = false)
    return _counttypes(t, input)
end

"""
    countsubcommunities(p::AbstractPartition)

Returns number of subcommunities in a partition, p.
"""
function countsubcommunities(p::AbstractPartition)
    return _countsubcommunities(p)
end

"""
    getnames(t::AbstractTypes)

Returns the names of the types of the AbstractTypes.
"""
function getnames(t::AbstractTypes, input::Bool = false)
    return _getnames(t, input)
end

"""
    getnames(p::AbstractPartition)

Returns the names of the subcommunities of the AbstractPartition.
"""
function getnames(p::AbstractPartition)
    return _getnames(p)
end

"""
    gettypenames(m::AbstractMetacommunity)

Returns the names of the types of the AbstractMetacommunity.
"""
function gettypenames(m::AbstractMetacommunity, input::Bool = false)
    return getnames(gettypes(m, input))
end

"""
    getsubcommunitynames(m::AbstractMetacommunity)

Returns the names of the subcommunities of the AbstractMetacommunity.
"""
function getsubcommunitynames(m::AbstractMetacommunity)
    return getnames(getpartition(m))
end

"""
    getabundance(m::AbstractMetacommunity)

Returns the abundances array of the metacommunity.
"""
function getabundance(m::AbstractMetacommunity, input::Bool = false)
    return _getabundance(m, input)
end

"""
    getmetaabundance(m::AbstractMetacommunity)

Returns the metacommunity abundances of the metacommunity.
"""
function getmetaabundance(m::AbstractMetacommunity, input::Bool = false)
    return _getmetaabundance(m, input)
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
    calcabundance(t::AbstractTypes, a::AbstractArray)

Calculates the abundance a for AbstractTypes, t (if necessary).
"""
function calcabundance(t::AbstractTypes, a::AbstractArray)
    return _calcabundance(t, a)
end

"""
    calcsimilarity(t::AbstractTypes, a::AbstractArray)

Retrieves (and possibly calculates) a similarity matrix from t.
"""
function calcsimilarity(t::AbstractTypes, a::AbstractArray)
    return _calcsimilarity(t, a)
end

"""
    calcordinariness(t::AbstractTypes, a::AbstractArray)

Calculates the ordinariness of abundance a from AbstractTypes, t.
"""
function calcordinariness(t::AbstractTypes, a::AbstractArray)
    return _calcordinariness(t, a)
end
