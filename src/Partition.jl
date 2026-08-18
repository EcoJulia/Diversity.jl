# SPDX-License-Identifier: BSD-2-Clause

"""
    Subcommunities(num)

AbstractPartition subtype with multiple subcommunities.

"""
struct Subcommunities <: Diversity.API.AbstractPartition{Nothing}
    num::Int64
    names::Vector{String}

    function Subcommunities(num::Integer)
        num > 0 || error("Too few subcommunities")
        return new(num, map(x -> "$x", 1:num))
    end

    function Subcommunities(names::Vector{String})
        num = length(names)
        num > 0 || error("Too few subcommunities")
        return new(num, names)
    end
end

import Diversity.API._getsubcommunitynames
function _getsubcommunitynames(sc::Subcommunities)
    return sc.names
end

import Diversity.API._countsubcommunities
function _countsubcommunities(sub::Subcommunities)
    return sub.num
end

"""
    Onecommunity

AbstractPartition subtype containing only one subcommunity.
"""
struct Onecommunity <: Diversity.API.AbstractPartition{Nothing}
    namev::Vector{String}

    function Onecommunity(name::String = "1")
        return new([name])
    end
end

function _getsubcommunitynames(oc::Onecommunity)
    return oc.namev
end

function _countsubcommunities(::Onecommunity)
    return 1
end

import Diversity.API._subsetpartition
# Selecting the one subcommunity there is leaves an undivided metacommunity undivided, rather than
# demoting it to a Subcommunities of one. The only other selection available from one subcommunity
# is none at all, which falls through to Subcommunities and is refused there -- a metacommunity with
# no subcommunities is not a metacommunity.
function _subsetpartition(oc::Onecommunity, idx)
    return length(idx) == 1 ? oc : Subcommunities(oc.namev[idx])
end
