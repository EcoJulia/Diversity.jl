using DataFrames

importall Diversity.API

"""
    UniqueTypes

A subtype of AbstractTypes where all individuals are completely
distinct. This type is the simplest AbstractTypes subtype, which
identifies all individuals as unique and completely distinct from each
other.

"""
immutable UniqueTypes <: AbstractTypes
    num::Int64
    names::Vector{String}

    function UniqueTypes(num::Integer)
        num > 0 || error("Too few species")
        new(num, map(x -> "$x", 1:num))
    end

    function UniqueTypes(names::Vector{String})
        num = length(names)
        num > 0 || error("Too few species")
        new(num, names)
    end
end

function _counttypes(ut::UniqueTypes, ::Bool)
    return ut.num
end

function _gettypenames(ut::UniqueTypes, ::Bool)
    return ut.names
end

function _calcsimilarity(ut::UniqueTypes, ::Real)
    return eye(ut.num)
end

function _calcordinariness(::UniqueTypes, abundances::AbstractArray, ::Real)
    return abundances
end

"""
    Species

A subtype of AbstractTypes where all species are completely distinct.
This type is the simplest AbstractTypes subtype, which identifies all
species as unique and completely distinct from each other.

"""
const Species = UniqueTypes

"""
    Taxonomy

A subtype of AbstractTypes with similarity between related taxa,
creating taxonomic similarity matrices.

"""
immutable Taxonomy{FP <: AbstractFloat} <: AbstractTypes
    speciesinfo::DataFrame
    taxa::Dict{Symbol, FP}
    typelabel::Symbol
    
    function (::Type{Taxonomy{FP}}){FP <: AbstractFloat}(speciesinfo::DataFrame,
                                                         taxa::Dict{Symbol, FP},
                                                         typelabel::Symbol)
        sort(speciesinfo.colindex.names) == sort([keys(taxa)...]) ||
            error("Taxon labels do not match similarity values")
        typelabel ∈ speciesinfo.colindex.names ||
            error("$typelabel not found in DataFrame column names")
        new{FP}(speciesinfo, taxa, typelabel)
    end
end

function Taxonomy(speciesinfo::DataFrame, taxa::Dict,
                  typelabel::Symbol = :Species)
    Taxonomy{valtype(taxa)}(speciesinfo, taxa, typelabel)
end

function floattypes{FP}(::Taxonomy{FP})
    return Set([FP])
end

function _counttypes(tax::Taxonomy, ::Bool)
    return nrow(tax.speciesinfo)
end

function _gettypenames(tax::Taxonomy, ::Bool)
    return tax.speciesinfo[tax.typelabel]
end

function _calcsimilarity(::Taxonomy, ::Real)
    error("Can't generate a taxonomic similarity matrix yet")
end

"""
    GeneralTypes{FP, M}

An AbstractTypes subtype with a general similarity matrix. This
subtype simply holds a matrix with similarities between individuals.

# Members:

- `z` A two-dimensional matrix representing similarity between
individuals.
"""
immutable GeneralTypes{FP <: AbstractFloat, M <: AbstractMatrix} <: AbstractTypes
    """
        z

    A two-dimensional matrix representing similarity between
    individuals.
    """
    z::M

    """
        names

    Optional vector of type names.
    """
    names::Vector{String}

    """
    # Constructor for GeneralTypes

    Creates an instance of the GeneralTypes class, with an arbitrary
    similarity matrix.
    """
    function (::Type{GeneralTypes{FP, M}}){FP <: AbstractFloat,
                                           M <: AbstractMatrix}(zmatrix::M)
        size(zmatrix, 1) == size(zmatrix, 2) ||
            throw(DimensionMismatch("Similarity matrix is not square"))

        minimum(zmatrix) ≥ 0 || throw(DomainError())
        maximum(zmatrix) ≤ 1 || warn("Similarity matrix has values above 1")

        new{FP, M}(zmatrix, map(x -> "$x", 1:size(zmatrix, 1)))
    end

    function (::Type{GeneralTypes{FP, M}}){FP <: AbstractFloat,
                                           M <: AbstractMatrix}(zmatrix::M,
                                                                names::Vector{String})
        size(zmatrix, 1) == size(zmatrix, 2) ||
            throw(DimensionMismatch("Similarity matrix is not square"))

        minimum(zmatrix) ≥ 0 || throw(DomainError())
        maximum(zmatrix) ≤ 1 || warn("Similarity matrix has values above 1")

        length(names) == size(zmatrix, 1) ||
            error("Species name vector does not match similarity matrix")

        new{FP, M}(zmatrix, names)
    end
end

function GeneralTypes{FP <: AbstractFloat}(zmatrix::AbstractMatrix{FP})
    GeneralTypes{FP, typeof(zmatrix)}(zmatrix)
end

function GeneralTypes{FP <: AbstractFloat}(zmatrix::AbstractMatrix{FP},
                                           names::Vector{String})
    GeneralTypes{FP, typeof(zmatrix)}(zmatrix, names)
end

function floattypes{FP, M}(::GeneralTypes{FP, M})
    return Set([FP])
end

function _counttypes(gt::GeneralTypes, ::Bool)
    return size(gt.z, 1)
end

function _gettypenames(gt::GeneralTypes, ::Bool)
    return gt.names
end

function _calcsimilarity(gt::GeneralTypes, ::Real)
    return gt.z
end
