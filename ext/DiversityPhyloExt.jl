module DiversityPhyloExt

isdefined(Base, :get_extension) ? (using Phylo) : (using ..Phylo)
import Diversity: AbstractPhyloTypes, PhyloBranches
using Diversity
using Diversity.API
using Statistics
using AxisArrays

import Diversity.API: _addedoutputcols
function _addedoutputcols(::Diversity.AbstractPhyloTypes{TS}) where
    {LABEL, RT, NL, N, B, TS <: TreeSet{LABEL, RT, NL, N, B, <: AbstractTree}}
    return Dict{Symbol, Type}(:treename => LABEL)
end

import Diversity.API: _getaddedoutput
function _getaddedoutput(pt::Diversity.AbstractPhyloTypes{TS}) where
    {LABEL, RT, NL, N, B, TS <: TreeSet{LABEL, RT, NL, N, B, <: AbstractTree}}
    return Dict{Symbol, LABEL}(:treename => first(gettreenames(pt.tree)))
end

function Diversity.PhyloBranches(tree::Tree) where Tree <: AbstractTree
    leafnames = [getnodename(tree, node)
                 for node in traversal(tree, preorder) if isleaf(tree, node)]
    nleaf = length(leafnames)
    nleaf > 0 || error("Too few species")
    leafinfo = Dict{String, Float64}()
    speciesinfo = Dict{String, Tuple{String, String, Float64}}()
    ancestralnames = String[]
    
    for leaf in leafnames
        branches = branchhistory(tree, leaf)
        leafinfo["$leaf"] = heighttoroot(tree, leaf)
        for branch in branches
            branchname = getbranchname(tree, branch)
            name = "$leaf : $branchname"
            push!(ancestralnames, name)
            speciesinfo[name] =
                tuple("$leaf", "$branchname",
                      getlength(tree, branch) / leafinfo["$leaf"])
        end
    end

    nancestral = length(ancestralnames)
    Lbar = Statistics.mean(collect(values(leafinfo)))
    ancestralmatrix = Matrix{Float64}(undef, nancestral, nleaf)
    fill!(ancestralmatrix, 0.0)    
    for i in 1:nancestral
        for j in 1:nleaf
            if speciesinfo[ancestralnames[i]][1] == leafnames[j]
                ancestralmatrix[i, j] =
                    speciesinfo[ancestralnames[i]][3] * leafinfo[leafnames[j]]
            end
        end
    end
    
    Zmatrix = Matrix{Float64}(undef, nancestral, nancestral)
    fill!(Zmatrix, 0.0)
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
    
    return Diversity.PhyloBranches{Tree}(tree, nleaf, nancestral, leafnames,
                                         ancestralnames, ancestralmatrix, Zmatrix)
end

function Diversity.PhyloBranches(treeset::TS) where TS <: TreeSet
    ntrees(treeset) == 1 ||
        error("Can currently only handle one tree in a PhyloSet")
    pt = Diversity.PhyloBranches(first(treeset))
    return Diversity.PhyloBranches{TS}(treeset, pt.nleaf, pt.nancestral,
                                       pt.leafnames, pt.ancestralnames,
                                       pt.ancestralmatrix, pt.Zmatrix)
end

import Diversity.API: _gettypenames
function _gettypenames(phy::Diversity.PhyloBranches, raw::Bool)
    return raw ? phy.leafnames : phy.ancestralnames
end

import Diversity.API: _counttypes
function _counttypes(phy::Diversity.PhyloBranches, raw::Bool)
    return raw ? phy.nleaf : phy.nancestral
end

import Diversity.API: _calcabundance
function _calcabundance(phy::Diversity.PhyloBranches,
                        raw::AbstractMatrix{<: AbstractFloat})
    processed = phy.ancestralmatrix * raw
    scale = sum(processed)
    processed ./= scale
    return processed, scale
end

import Diversity.API: _calcsimilarity
function _calcsimilarity(phy::Diversity.PhyloBranches, scale::Real)
    return phy.Zmatrix .* scale
end

import Diversity.API: _calcordinariness
function _calcordinariness(phy::Diversity.PhyloBranches,
                           processed::AbstractMatrix{<: AbstractFloat},
                           scale::Real)
    return (phy.Zmatrix * processed) .* scale
end

import Diversity.API: _getdiversityname
_getdiversityname(::Diversity.PhyloBranches) = "Phylogenetic Branch"

end
