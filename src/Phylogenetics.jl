using PhyloTrees
using Diversity
import Diversity.counttypes, Diversity.getsimilarity
import Diversity.getordinariness, Diversity.getnames

type Phylogeny{ND} <: AbstractTypes
    num::Int64
    names::Vector{String}
    phylo::NodeTree{ND}

    function (::Type{Phylogeny{ND}}){ND}(phylo::NodeTree{ND})
        names = getleafnames(phylo)
        num = length(names)
        num > 0 || error("Too few species")
        new{ND}(num, names, phylo)
    end
end

Phylogeny{ND}(phylo::NodeTree{ND}) = Phylogeny{ND}(phylo)
              
function counttypes(phy::Phylogeny)
    return phy.num
end

function getsimilarity(phy::Phylogeny)
    return eye(phy.num)
end

function getordinariness(::Phylogeny, abundances::AbstractArray)
    return abundances
end

function getnames(phy::Phylogeny)
    return phy.names
end
