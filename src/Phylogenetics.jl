using PhyloTrees
using Diversity
import Diversity.counttypes, Diversity.getsimilarity, Diversity.getabundance
import Diversity.getordinariness, Diversity.getnames

type Phylogeny{Tree <: AbstractTree} <: AbstractTypes
    tree::Tree
    nleaf::Int64
    nancestral::Int64
    leafnames::Vector{String}
    speciesnames::Vector{String}
    ancestralmatrix::Matrix{Float64}
    Zmatrix::Matrix{Float64}

    function (::Type{Phylogeny{Tree}}){Tree <: AbstractTree}(tree::Tree)
        leafnames = getleafnames(tree)
        nleaf = length(leafnames)
        nleaf > 0 || error("Too few species")
        leafinfo = Dict{String, Float64}()
        speciesinfo = Dict{String, Tuple{String, String, Float64}}()
        for leaf in leafnames
            branches = branchpath(tree, leaf)
            leafinfo["$leaf"] = getrootdistance(tree, leaf)
            for branch in branches
                name = "$leaf : $branch"
                speciesinfo[name] = tuple("$leaf", "$branch",
                                          getlength(tree, branch) / leafinfo["$leaf"])
            end
        end
        speciesnames = collect(keys(speciesinfo))
        nancestral = length(speciesnames)
        Lbar = mean(values(leafinfo))
        ancestralmatrix = Matrix{Float64}(nancestral, nleaf)
        fill!(ancestralmatrix, 0)
        for i in 1:nancestral
            for j in 1:nleaf
                if speciesinfo[speciesnames[i]][1] == leafnames[j]
                    ancestralmatrix[i, j] =
                        speciesinfo[speciesnames[i]][3] * leafinfo[leafnames[j]]
                end
            end
        end
        Zmatrix = Matrix{Float64}(nancestral, nancestral)
        fill!(Zmatrix, 0)
        for i in 1:nancestral
            for j in 1:nancestral
                if haskey(speciesinfo,
                          "$(speciesinfo[speciesnames[j]][1]) : $(speciesinfo[speciesnames[i]][2])")
                    Zmatrix[i, j] = 1.0 / leafinfo[speciesinfo[speciesnames[j]][1]]
                end
            end
        end
        new{Tree}(tree, nleaf, nancestral, leafnames, speciesnames,
                  ancestralmatrix, Zmatrix)
    end
end

Phylogeny{Tree <: AbstractTree}(tree::Tree) = Phylogeny{Tree}(tree)
              
function counttypes(phy::Phylogeny)
    return phy.nancestral
end

function getsimilarity(phy::Phylogeny)
    return phy.Zmatrix * phy.ancestralmatrix
end

function getordinariness(phy::Phylogeny, abundances::AbstractArray)
    return phy.Zmatrix * phy.ancestralmatrix * abundances
end

function getabundance(phy::Phylogeny, abundances::AbstractArray)
    return phy.ancestralmatrix * abundances
end

function getnames(phy::Phylogeny)
    return phy.speciesnames
end
