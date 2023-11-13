using DataFrames
using LinearAlgebra

"""
    UniqueTypes

A subtype of AbstractTypes where all individuals are completely
distinct. This type is the simplest AbstractTypes subtype, which
identifies all individuals as unique and completely distinct from each
other.

"""
struct UniqueTypes <: Diversity.API.AbstractTypes
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

import Diversity.API._hassimilarity
_hassimilarity(::UniqueTypes) = false

import Diversity.API._counttypes
function _counttypes(ut::UniqueTypes, ::Bool)
    return ut.num
end

import Diversity.API._gettypenames
function _gettypenames(ut::UniqueTypes, ::Bool)
    return ut.names
end

import Diversity.API._calcsimilarity
function _calcsimilarity(ut::UniqueTypes, ::Real)
    return I
end

import Diversity.API._calcordinariness
function _calcordinariness(::UniqueTypes, abundances::AbstractArray, ::Real)
    return abundances
end

import Diversity.API._getdiversityname
_getdiversityname(::UniqueTypes) = "Unique"

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
struct Taxonomy{FP <: AbstractFloat} <: Diversity.API.AbstractTypes
    speciesinfo::DataFrame
    taxa::Dict{Symbol, FP}
    typelabel::Symbol

    function Taxonomy{FP}(speciesinfo::DataFrame,
                          taxa::Dict{Symbol, FP},
                          typelabel::Symbol) where FP <: AbstractFloat
        sort(describe(speciesinfo)[!,:variable]) == sort([keys(taxa)...]) ||
            error("Taxon labels do not match similarity values")
        typelabel ∈ describe(speciesinfo)[!,:variable] ||
            error("$typelabel not found in DataFrame column names")
        new{FP}(speciesinfo, taxa, typelabel)
    end
end

function Taxonomy(speciesinfo::DataFrame, taxa::Dict,
                  typelabel::Symbol = :Species)
    Taxonomy{valtype(taxa)}(speciesinfo, taxa, typelabel)
end

import Diversity.API.floattypes
function floattypes(::Taxonomy{FP}) where FP <: AbstractFloat
    return Set([FP])
end

function _counttypes(tax::Taxonomy, ::Bool)
    return nrow(tax.speciesinfo)
end

function _gettypenames(tax::Taxonomy, ::Bool)
    return tax.speciesinfo[!,tax.typelabel]
end

function _calcsimilarity(::Taxonomy, ::Real)
    error("Can't generate a taxonomic similarity matrix yet")
end

_getdiversityname(::Taxonomy) = "Taxonomy"

"""
    GeneralTypes{FP, M, LABELS}

An AbstractTypes subtype with a general similarity matrix. This
subtype simply holds a matrix with similarities between individuals.

# Members:

- `z` A two-dimensional matrix representing similarity between
individuals.

- `names` Vector of type names.

"""
struct GeneralTypes{FP <: AbstractFloat,
                    M <: AbstractMatrix{FP},
                    LABELS <: AbstractVector} <: Diversity.API.AbstractTypes
    z::M
    names::LABELS

    function GeneralTypes(zmatrix::M, names::LABELS) where
        {FP <: AbstractFloat, M <: AbstractMatrix{FP}, LABELS <: AbstractVector}

        size(zmatrix, 1) == size(zmatrix, 2) ||
            throw(DimensionMismatch("Similarity matrix is not square"))

        minimum(zmatrix) ≥ 0 || throw(DomainError(minimum(zmatrix),
                                      "Similarities must be ≥ 0"))
        maximum(zmatrix) ≤ 1 || @warn "Similarity matrix has values above 1"

        length(names) == size(zmatrix, 1) ||
            error("Species name vector does not match similarity matrix")

        return new{FP, M, LABELS}(zmatrix, names)
    end    
end

"""
    GeneralTypes(zmatrix::M)
    GeneralTypes(zmatrix::M, names::LABELS)

Constructors for GeneralTypes. Creates an instance of the GeneralTypes class,
with an arbitrary `zmatrix` similarity matrix and an optional vector of type `names`.
"""
function GeneralTypes(zmatrix::M) where {FP <: AbstractFloat,
                                         M <: AbstractMatrix{FP}}
    ax = collect(axes(zmatrix, 1))
    return GeneralTypes(zmatrix, ax)
end

function floattypes(::GeneralTypes{FP}) where {FP <: AbstractFloat}
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

_getdiversityname(::GeneralTypes) = "Arbitrary Z"
