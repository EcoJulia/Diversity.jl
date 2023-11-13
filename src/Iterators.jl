using Diversity.API

import Base.IteratorSize, Base.length, Base.IteratorEltype, Base.HasEltype, Base.HasLength, Base.eltype

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

import Base.iterate
function iterate(ti::TypeIterator, state = 1)
    if state > counttypes(ti.metacommunity)
        return nothing
    end
    return ti.viewfn(ti.metacommunity, state), state + 1
end

function IteratorSize(::Type{<:TypeIterator})
    return HasLength()
end

function length(ti::TypeIterator)
    return counttypes(ti.metacommunity)
end

function IteratorEltype(::Type{<:TypeIterator})
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

function iterate(si::SubcommunityIterator, state = 1)
    if state > countsubcommunities(si.metacommunity)
        return nothing
    end
    return si.viewfn(si.metacommunity, state), state + 1
end

function IteratorSize(::Type{<:SubcommunityIterator})
    return HasLength()
end

function length(si::SubcommunityIterator)
    return countsubcommunities(si.metacommunity)
end

function IteratorEltype(::Type{<:SubcommunityIterator})
    return HasEltype()
end

function eltype(si::SubcommunityIterator)
    return eltype(si.viewfn(si.metacommunity, 1))
end
