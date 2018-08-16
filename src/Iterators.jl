using Compat
using Diversity.API

import Compat.IteratorSize, Base.length, Compat.IteratorEltype, Base.eltype

abstract type AbstractIterator{M <: AbstractMetacommunity} end

struct TypeIterator{M} <: AbstractIterator{M}
    metacommunity::M
    viewfn::Function
end

function TypeIterator(fn::Function, meta::M) where M <: AbstractMetacommunity
    n = ndims(fn(meta))
    if n == 2
        viewfn = (meta, i) -> view(fn(meta), i, :)
    elseif n == 1
        viewfn = (meta, i) -> view(fn(meta), i)
    else
        error("Can't iterate over types for the function of the metacommunity - $fn")
    end
    return TypeIterator(meta, viewfn)
end

function TypeIterator(meta::M) where M <: AbstractMetacommunity
    return TypeIterator(getabundance, meta)
end

if VERSION < v"0.7.0-"
import Base.start, Base.next, Base.done
function start(ti::TypeIterator)
    return 1
end

function next(ti::TypeIterator, state)
    return (ti.viewfn(ti.metacommunity, state), state + 1)
end

function done(ti::TypeIterator, state)
    return state > counttypes(ti.metacommunity)
end

else
import Base.iterate
function iterate(ti::TypeIterator, state = 1)
    if state > counttypes(ti.metacommunity)
        return nothing
    end
    return ti.viewfn(ti.metacommunity, state), state + 1
end

end

function IteratorSize(ti::Type{TypeIterator})
    return HasLength()
end

function length(ti::TypeIterator)
    return counttypes(ti.metacommunity)
end

function IteratorEltype(ti::Type{TypeIterator})
    return HasEltype()
end

function eltype(ti::TypeIterator)
    return eltype(ti.viewfn(ti.metacommunity, 1))
end

struct SubcommunityIterator{M} <: AbstractIterator{M}
    metacommunity::M
    viewfn::Function
end

function SubcommunityIterator(fn::Function,
                              meta::M) where M <: AbstractMetacommunity
    n = ndims(fn(meta))
    if n == 2
        viewfn = (meta, i) -> view(fn(meta), :, i)
    elseif n == 1
        viewfn = (meta, i) -> view(fn(meta), i)
    else
        error("Can't iterate over subcommunities for the " *
              "function of the metacommunity - $fn")
    end
    return SubcommunityIterator(meta, viewfn)
end

function SubcommunityIterator(meta::M) where M <: AbstractMetacommunity
    return SubcommunityIterator(getabundance, meta)
end

if VERSION < v"0.7.0-"
function start(si::SubcommunityIterator)
    return 1
end

function next(si::SubcommunityIterator, state)
    return (si.viewfn(si.metacommunity, state), state + 1)
end

function done(si::SubcommunityIterator, state)
    return state > countsubcommunities(si.metacommunity)
end

else

function iterate(si::SubcommunityIterator, state = 1)
    if state > countsubcommunities(si.metacommunity)
        return nothing
    end
    return si.viewfn(si.metacommunity, state), state + 1
end

end


function IteratorSize(si::Type{SubcommunityIterator})
    return HasLength()
end

function length(si::SubcommunityIterator)
    return countsubcommunities(si.metacommunity)
end

function IteratorEltype(si::Type{SubcommunityIterator})
    return HasEltype()
end

function eltype(si::SubcommunityIterator)
    return eltype(si.viewfn(si.metacommunity, 1))
end
