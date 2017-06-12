using Phylo
importall Diversity.API

type PhyloTypes{Tree <: AbstractTree} <: AbstractTypes
    tree::Tree
    nleaf::Int64
    nancestral::Int64
    leafnames::Vector{String}
    ancestralnames::Vector{String}
    ancestralmatrix::Matrix{Float64}
    Zmatrix::Matrix{Float64}

    function (::Type{PhyloTypes{Tree}}){Tree <: AbstractTree}(tree::Tree)
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
                          "$(speciesinfo[ancestralnames[j]][1]) : $(speciesinfo[ancestralnames[i]][2])")
                    Zmatrix[i, j] = 1.0 / leafinfo[speciesinfo[ancestralnames[j]][1]]
                end
            end
        end
        new{Tree}(tree, nleaf, nancestral, leafnames, ancestralnames,
                  ancestralmatrix, Zmatrix)
    end
end

PhyloTypes{Tree <: AbstractTree}(tree::Tree) = PhyloTypes{Tree}(tree)

function _getnames(phy::PhyloTypes, input::Bool)
    return input ? phy.leafnames : phy.ancestralnames
end

function _counttypes(phy::PhyloTypes, input::Bool)
    return input ? phy.nleaf : phy.nancestral
end

function _calcabundance(phy::PhyloTypes, abundances::AbstractArray)
    asabundances = phy.ancestralmatrix * abundances
    asabundances /= sum(asabundances)
    return asabundances
end

function _calcsimilarity(phy::PhyloTypes, abundances::AbstractArray)
    asabundances = phy.ancestralmatrix * abundances
    return phy.Zmatrix * sum(asabundances)
end

function _calcordinariness(phy::PhyloTypes, abundances::AbstractArray)
    return phy.Zmatrix * phy.ancestralmatrix * abundances
end
