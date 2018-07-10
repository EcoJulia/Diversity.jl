using Phylo
using Diversity

import Diversity.Phylogenetics: PhyloTypes
struct PhyloTypes{Tree <: AbstractTree} <: Diversity.API.AbstractTypes
    tree::Tree
    nleaf::Int64
    nancestral::Int64
    leafnames::Vector{String}
    ancestralnames::Vector{String}
    ancestralmatrix::Matrix{Float64}
    Zmatrix::Matrix{Float64}

    function PhyloTypes(tree::Tree) where Tree <: AbstractTree
        leafnames = getleafnames(tree)
        nleaf = length(leafnames)
        nleaf > 0 || error("Too few species")
        leafinfo = Dict{String, Float64}()
        speciesinfo = Dict{String, Tuple{String, String, Float64}}()
        for leaf in leafnames
            branches = branchhistory(tree, leaf)
            leafinfo["$leaf"] = heighttoroot(tree, leaf)
            for branch in branches
                name = "$leaf : $branch"
                speciesinfo[name] = tuple("$leaf", "$branch",
                                          getlength(tree, branch) / leafinfo["$leaf"])
            end
        end
        ancestralnames = collect(keys(speciesinfo))
        nancestral = length(ancestralnames)
        Lbar = mean(values(leafinfo))
        ancestralmatrix = Matrix{Float64}(nancestral, nleaf)
        fill!(ancestralmatrix, 0)
        for i in 1:nancestral
            for j in 1:nleaf
                if speciesinfo[ancestralnames[i]][1] == leafnames[j]
                    ancestralmatrix[i, j] =
                        speciesinfo[ancestralnames[i]][3] * leafinfo[leafnames[j]]
                end
            end
        end
        Zmatrix = Matrix{Float64}(nancestral, nancestral)
        fill!(Zmatrix, 0)
        for i in 1:nancestral
            for j in 1:nancestral
                if haskey(speciesinfo,
                          "$(speciesinfo[ancestralnames[j]][1]) : " *
                          "$(speciesinfo[ancestralnames[i]][2])")
                    Zmatrix[i, j] = 1.0 /
                        leafinfo[speciesinfo[ancestralnames[j]][1]]
                end
            end
        end
        new{Tree}(tree, nleaf, nancestral, leafnames, ancestralnames,
                  ancestralmatrix, Zmatrix)
    end
end

import Diversity.API._gettypenames
function _gettypenames(phy::PhyloTypes, raw::Bool)
    return raw ? phy.leafnames : phy.ancestralnames
end

import Diversity.API._counttypes
function _counttypes(phy::PhyloTypes, raw::Bool)
    return raw ? phy.nleaf : phy.nancestral
end

import Diversity.API._calcabundance
function _calcabundance(phy::PhyloTypes, raw::M) where {FP <: AbstractFloat,
                                                        M <: AbstractMatrix{FP}}
    processed = phy.ancestralmatrix * raw
    scale = sum(processed)
    processed ./= scale
    return processed, scale
end

import Diversity.API._calcsimilarity
function _calcsimilarity(phy::PhyloTypes, scale::R) where R <: Real
    return phy.Zmatrix .* scale
end

import Diversity.API._calcordinariness
function _calcordinariness(phy::PhyloTypes, processed::M, scale::R) where
    {FP <: AbstractFloat, M <: AbstractMatrix{FP}, R <: Real}
    return (phy.Zmatrix * processed) .* scale
end
