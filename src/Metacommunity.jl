using Compat
using DataFrames

"""
    Subcommunities(numsub)

AbstractPartition subtype with multiple subcommunities.

"""
immutable Subcommunities <: AbstractPartition
    numsub::Int64

    function Subcommunities(numsub::Integer)
        numsub > 0 || error("Too few subcommunities")
        new(numsub)
    end
end

function countsubcommunities(sub::Subcommunities)
    return sub.numsub
end

"""
    Onecommunity

AbstractPartition subtype containing only one subcommunity.
"""
immutable Onecommunity <: AbstractPartition end

function countsubcommunities(::Onecommunity)
    return 1
end

"""
    UniqueTypes

A subtype of AbstractTypes where all individuals are completely
distinct. This type is the simplest AbstractTypes subtype, which
identifies all individuals as unique and completely distinct from each
other.

"""
immutable UniqueTypes <: AbstractTypes
    num::Int64
    names::Nullable{Vector{String}}

    function UniqueTypes(num::Integer)
        num > 0 || error("Too few species")
        new(num, Nullable{Vector{String}}())
    end

    function UniqueTypes(names::Vector{String})
        num = length(names)
        num > 0 || error("Too few species")
        new(num, Nullable(names))
    end
end

function counttypes(ut::UniqueTypes)
    return ut.num
end

function getsimilarity(ut::UniqueTypes)
    return eye(ut.num)
end

function getordinariness(::UniqueTypes, abundances::AbstractArray)
    return abundances
end

function hasnames(ut::UniqueTypes)
    return !isnull(ut.names)
end
                  
function getnames(ut::UniqueTypes)
    return get(ut.names)
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
    
    function (::Type{Taxonomy{FP}}){FP <:
        AbstractFloat}(speciesinfo::DataFrame,
                       taxa::Dict{Symbol, FP})
        sort(speciesinfo.colindex.names) == sort([keys(taxa)...]) ||
        error("Taxon labels do not match similarity values")
        new{FP}(speciesinfo, taxa)
    end
end

function Taxonomy(speciesinfo::DataFrame, taxa::Dict)
    Taxonomy{valtype(taxa)}(speciesinfo, taxa)
end

function counttypes(tax::Taxonomy)
    return nrow(tax.speciesinfo)
end

function getsimilarity(::Taxonomy)
    error("Can't generate a taxonomic similarity matrix yet")
end

function floattypes{FP}(::Taxonomy{FP})
    return Set([FP])
end

function hasnames(tax::Taxonomy)
    return :Species ∈ tax.speciesinfo.colindex.names
end
                  
function getnames(tax::Taxonomy)
    return tax.speciesinfo[:Species]
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
    names::Nullable{Vector{String}}

    """
    # Constructor for GeneralTypes

    Creates an instance of the GeneralTypes class, with an arbitrary similarity matrix.

    # Arguments:
    - `zmatrix`: similarity matrix
    """
    function (::Type{GeneralTypes{FP, M}}){FP <: AbstractFloat,
        M <: AbstractMatrix}(zmatrix::M)
        size(zmatrix, 1) == size(zmatrix, 2) ||
        throw(DimensionMismatch("Similarity matrix is not square"))

        minimum(zmatrix) ≥ 0 || throw(DomainError())
        maximum(zmatrix) ≤ 1 || warn("Similarity matrix has values above 1")

        new{FP, M}(zmatrix, Nullable{Vector{String}}())
    end

    function (::Type{GeneralTypes{FP, M}}){FP <: AbstractFloat,
        M <: AbstractMatrix}(zmatrix::M, names::Vector{String})
        size(zmatrix, 1) == size(zmatrix, 2) ||
        throw(DimensionMismatch("Similarity matrix is not square"))

        minimum(zmatrix) ≥ 0 || throw(DomainError())
        maximum(zmatrix) ≤ 1 || warn("Similarity matrix has values above 1")

        length(names) == size(zmatrix, 1) ||
        error("Species name vector does not match similarity matrix")

        new{FP, M}(zmatrix, Nullable(names))
    end
end

function GeneralTypes{FP <: AbstractFloat}(zmatrix::AbstractMatrix{FP})
    GeneralTypes{FP, typeof(zmatrix)}(zmatrix)
end

function GeneralTypes{FP <: AbstractFloat}(zmatrix::AbstractMatrix{FP}, names::Vector{String})
    GeneralTypes{FP, typeof(zmatrix)}(zmatrix, names)
end

function counttypes(gt::GeneralTypes)
    return size(gt.z, 1)
end

function getsimilarity(gt::GeneralTypes)
    return gt.z
end

function floattypes{FP, M}(::GeneralTypes{FP, M})
    return Set([FP])
end

function hasnames(gt::GeneralTypes)
    return !isnull(gt.names)
end
                  
function getnames(gt::GeneralTypes)
    return get(gt.names)
end

"""
    Metacommunity{FP, A, Part, Sim}

Metacommunity type, representing a whole metacommunity containing a
single community or a collection of subcommunities. The metacommunity
of individuals *may* be further partitioned into smaller groups. For
instance this may be an ecosystem, which consists of a series of
subcommunities. The AbstractPartition subtype within it stores
relative abundances of different types, e.g. species, and also allows
for similarity between individuals.

# Constructor:

**Metacommunity(abundances::AbstractArray, part::AbstractPartition, types::AbstractTypes)**

# Members:

- `abundances` the abundance matrix for the metacommunity.

- `partition` the instance of the AbstractPartition subtype, containing the
  subcommunities.

- `types` The instance of the AbstractTypes subtype, from which
  similarities between individuals can be calculated.

- `ordinariness` A cache of the ordinariness of the individuals in the
  Partition. Should only be accessed through
  getordinariness!(::Metacommunity), which will populate the cache if
  it has not yet been calculated.

"""
type Metacommunity{FP, A, Sim, Part} <: AbstractMetacommunity{FP, A, Sim, Part}
    abundances::A
    types::Sim
    partition::Part
    ordinariness::Nullable{A}

    function (::Type{Metacommunity{FP, A, Sim, Part}}){FP <: AbstractFloat,
        A <: AbstractArray,
        Sim <: AbstractTypes,
        Part <: AbstractPartition}(abundances::A, types::Sim, part::Part)
        mcmatch(abundances, types, part) ||
        throw(ErrorException("Type or size mismatch between abundance array, partition and type list"))
        new{FP, A, Sim, Part}(abundances, types, part, Nullable{A}())
    end

    function (::Type{Metacommunity{FP, A, Sim, Part}}){FP <: AbstractFloat,
        A <: AbstractArray,
        Sim <: AbstractTypes,
        Part <: AbstractPartition}(abundances::A, meta::Metacommunity{FP, A, Sim, Part})
        mcmatch(abundances, meta.types, meta.part) ||
        throw(ErrorException("Type or size mismatch between abundance array, partition and type list"))
        new{FP, A, Sim, Part}(abundances, meta.types, meta.part, meta.ordinariness)
    end
end

function Metacommunity{A <: AbstractArray,
    Meta <: AbstractMetacommunity}(abundances::A, meta::Meta)
    types = gettypes(meta)
    part = getpartition(meta)
    Metacommunity{eltype{A}, A, typeof(types), typeof(part)}(abundances, types, part)
end

function Metacommunity{M <: AbstractMatrix, Sim <: AbstractTypes,
    Part <: AbstractPartition}(abundances::M,
                               types::Sim = UniqueTypes(size(abundances, 1)),
                               part::Part = Subcommunities(size(abundances, 2)))
    Metacommunity{eltype(M), M, Sim, Part}(abundances, types, part)
end

function Metacommunity{V <: AbstractVector, Sim <: AbstractTypes,
    Part <: AbstractPartition}(abundances::V,
                               types::Sim = UniqueTypes(size(abundances, 1)),
                               part::Part = Onecommunity())
    Metacommunity{eltype(V), V, Sim, Part}(abundances, types, part)
end

function Metacommunity{V <: AbstractVector, M <: AbstractMatrix}(abundances::V, zmatrix::M)
    Metacommunity{eltype(V), V,
    GeneralTypes, Onecommunity}(abundances, GeneralTypes(zmatrix), Onecommunity())
end

function Metacommunity{M <: AbstractMatrix}(abundances::M, zmatrix::M)
    Metacommunity{eltype(M), M,
    GeneralTypes, Subcommunities}(abundances, GeneralTypes(zmatrix),
                                  Subcommunities(size(abundances, 2)))
end

function getabundance{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.abundances
end

function gettypes{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.types
end

function getpartition{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.partition
end

@inline function getordinariness!(meta::Metacommunity)
    if isnull(meta.ordinariness)
        meta.ordinariness = getordinariness(meta.types, meta.abundances)
    end
    get(meta.ordinariness)
end
