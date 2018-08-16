using .Phylo
using Diversity
using Diversity.API
using Compat.Statistics

struct PhyloTypes{Tree <: AbstractTree} <: Diversity.API.AbstractTypes
    tree::Tree
    nleaf::Int64
    nancestral::Int64
    leafnames::Vector{String}
    ancestralnames::Vector{String}
    ancestralmatrix::Matrix{Float64}
    Zmatrix::Matrix{Float64}
end

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
    Lbar = Compat.Statistics.mean(collect(values(leafinfo)))
    ancestralmatrix = fill(zero(Float64), (nancestral, nleaf))
    for i in 1:nancestral
        for j in 1:nleaf
            if speciesinfo[ancestralnames[i]][1] == leafnames[j]
                ancestralmatrix[i, j] =
                    speciesinfo[ancestralnames[i]][3] * leafinfo[leafnames[j]]
            end
        end
    end
    Zmatrix = fill(zero(Float64), (nancestral, nancestral))
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
    return PhyloTypes{Tree}(tree, nleaf, nancestral,
                            leafnames, ancestralnames,
                            ancestralmatrix, Zmatrix)
end

function PhyloTypes(treeset::TS) where TS <: TreeSet
    ntrees(treeset) == 1 ||
        error("Can currently only handle one tree in a PhyloSet")
    pt = PhyloTypes(first(treeset))
    return PhyloTypes{TS}(treeset, pt.nleaf, pt.nancestral,
                          pt.leafnames, pt.ancestralnames,
                          pt.ancestralmatrix, pt.Zmatrix)
end

import Diversity.API: _gettypenames
function _gettypenames(phy::Diversity.PhyloTypes, raw::Bool)
    return raw ? phy.leafnames : phy.ancestralnames
end

import Diversity.API: _counttypes
function _counttypes(phy::Diversity.PhyloTypes, raw::Bool)
    return raw ? phy.nleaf : phy.nancestral
end

import Diversity.API: _calcabundance
function _calcabundance(phy::Diversity.PhyloTypes, raw::M) where
    {M <: AbstractMatrix{<:AbstractFloat}}
    processed = phy.ancestralmatrix * raw
    scale = sum(processed)
    processed ./= scale
    return processed, scale
end

import Diversity.API: _calcsimilarity
function _calcsimilarity(phy::Diversity.PhyloTypes, scale::R) where R <: Real
    return phy.Zmatrix .* scale
end

import Diversity.API: _calcordinariness
function _calcordinariness(phy::Diversity.PhyloTypes, processed::M, scale::R) where
    {FP <: AbstractFloat, M <: AbstractMatrix{FP}, R <: Real}
    return (phy.Zmatrix * processed) .* scale
end

import Diversity.API: _getdiversityname
_getdiversityname(::Diversity.PhyloTypes) = "Phylogenetic"

import Diversity.API: _addedoutputcols
_addedoutputcols(::Diversity.PhyloTypes{TS}) where
    {LABEL, NL, BL, TS <: TreeSet{LABEL, NL, BL, <: AbstractTree}} =
    Dict{Symbol, Type}(:treename => LABEL)

import Diversity.API: _getaddedoutput
_getaddedoutput(pt::Diversity.PhyloTypes{TS}) where
    {LABEL, NL, BL, TS <: TreeSet{LABEL, NL, BL, <: AbstractTree}} =
    Dict{Symbol, LABEL}(:treename => first(treenameiter(pt.tree)))
