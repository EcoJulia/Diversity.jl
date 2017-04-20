using PhyloTrees
using Diversity
import Diversity.counttypes, Diversity.getsimilarity
import Diversity.getordinariness, Diversity.getnames

type Phylogeny{Tree <: AbstractTree} <: AbstractTypes
    num::Int64
    leafnames::Vector{String}
    tree::Tree
    speciesnames::Vector{String}
    ancestralmatrix::Matrix{Float64}
    Zmatrix::Matrix{Float64}

    function (::Type{Phylogeny{Tree}}){Tree <: AbstractTree}(tree::Tree)
        leafnames = getleafnames(tree)
        num = length(leafnames)
        num > 0 || error("Too few species")
        leafinfo = Dict{String, Float64}()
        speciesinfo = Dict{String, Tuple{String, Float64}}()
        for leaf in leafnames
            branches = branchpath(tree, leaf)
            leafinfo[leaf] = getrootdistance(tree, leaf)
            for branch in branches
                name = "$leaf : $branch"
                speciesinfo[name] = tuple(leaf, getlength(tree, branch))
            end
        end
        
        new{Tree}(num, leafnames, tree, collect(keys(speciesinfo)),
                  Matrix{Float64}(), Matrix{Float64}())
    end
end

Phylogeny{Tree <: AbstractTree}(tree::Tree) = Phylogeny{Tree}(tree)
              
function counttypes(phy::Phylogeny)
    return phy.num
end

function getsimilarity(phy::Phylogeny)
    return eye(phy.num)
end

function getordinariness(::Phylogeny, abundances::AbstractArray)
    return Zmatrix * abundances
end

function getabundance(::Phylogeny, abundances::AbstractArray)
    return ancestralmatrix * abundances
end

function getnames(phy::Phylogeny)
    return phy.leafnames
end
