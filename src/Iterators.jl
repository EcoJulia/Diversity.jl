using Diversity.API

import Base.start, Base.next, Base.done
import Base.iteratorsize, Base.length, Base.iteratoreltype, Base.eltype

# length of an AbstractTypes is the number of types
function length(t::AbstractTypes)
    return _counttypes(t)
end

# length of an AbstractPartition is the number of subcommunities
function length(p::AbstractPartition)
    return _countsubcommunities(p)
end

immutable TypeIterator{M <: AbstractMetacommunity}
    metacommunity::M
    viewfn::Function

    (::Type{TypeIterator{M}}){M}(metacommunity::M,
                                 viewfn::Function) = new{M}(metacommunity, viewfn)    
end

function TypeIterator(meta::AbstractMetacommunity,
                      fn::Function = getabundance)
    n = ndims(fn(meta))
    if n == 2
        viewfn = (meta, i) -> view(fn(meta), i, :)
    elseif n == 1
        viewfn = (meta, i) -> view(fn(meta), i)
    else
        error("Can't iterate over types for the function of the metacommunity - $fn")
    end
    return TypeIterator{typeof(meta)}(meta, viewfn)
end
 
function start(ti::TypeIterator)
    return 1
end

function next(ti::TypeIterator, state)
    return (ti.viewfn(ti.metacommunity, state), state + 1)
end

function done(ti::TypeIterator, state)
    return state > counttypes(ti.metacommunity)
end

function iteratorsize(ti::Type{TypeIterator})
    return HasLength()
end

function length(ti::TypeIterator)
    return counttypes(ti.metacommunity)
end

function iteratoreltype(ti::Type{TypeIterator})
    return HasEltype()
end

function eltype{T <: AbstractMetacommunity}(ti::TypeIterator{T})
    return eltype(ti.viewfn(ti.metacommunity, 1))
end



immutable SubcommunityIterator{M <: AbstractMetacommunity}
    metacommunity::M
    viewfn::Function

    (::Type{SubcommunityIterator{M}}){M}(metacommunity::M,
                                         viewfn::Function) = new{M}(metacommunity, viewfn)
end

function SubcommunityIterator(meta::AbstractMetacommunity,
                              fn::Function = getabundance)
    n = ndims(fn(meta))
    if n == 2
        viewfn = (meta, i) -> view(fn(meta), :, i)
    elseif n == 1
        viewfn = (meta, i) -> view(fn(meta), i)
    else
        error("Can't iterate over subcommunities for the function of the metacommunity - $fn")
    end
    return SubcommunityIterator{typeof(meta)}(meta, viewfn)
end
 
function start(si::SubcommunityIterator)
    return 1
end

function next(si::SubcommunityIterator, state)
    return (si.viewfn(si.metacommunity, state), state + 1)
end

function done(si::SubcommunityIterator, state)
    return state > countsubcommunities(si.metacommunity)
end

function iteratorsize(si::Type{SubcommunityIterator})
    return HasLength()
end

function length(si::SubcommunityIterator)
    return countsubcommunities(si.metacommunity)
end

function iteratoreltype(si::Type{SubcommunityIterator})
    return HasEltype()
end

function eltype{T <: AbstractMetacommunity}(si::SubcommunityIterator{T})
    return eltype(si.viewfn(si.metacommunity, 1))
end
