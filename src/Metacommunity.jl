using Compat
using DataFrames

importall Diversity.API

"""
    Subcommunities(num)

AbstractPartition subtype with multiple subcommunities.

"""
immutable Subcommunities <: AbstractPartition
    num::Int64
    names::Vector{String}

    function Subcommunities(num::Integer)
        num > 0 || error("Too few subcommunities")
        new(num, map(x -> "$x", 1:num))
    end

    function Subcommunities(names::Vector{String})
        num = length(names)
        num > 0 || error("Too few subcommunities")
        new(num, names)
    end
end

function _countsubcommunities(sub::Subcommunities)
    return sub.num
end

function _getnames(sc::Subcommunities)
    return sc.names
end

"""
    Onecommunity

AbstractPartition subtype containing only one subcommunity.
"""
immutable Onecommunity <: AbstractPartition
    name::Vector{String}

    function Onecommunity(name::String = "1")
        new([name])
    end
end

function _countsubcommunities(::Onecommunity)
    return 1
end

function _getnames(oc::Onecommunity)
    return oc.name
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

function _getnames(ut::UniqueTypes, ::Bool)
    return ut.names
end
    
function _calcsimilarity(ut::UniqueTypes, ::AbstractArray)
    return eye(ut.num)
end

function _calcordinariness(::UniqueTypes, abundances::AbstractArray)
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

function _getnames(tax::Taxonomy, ::Bool)
    return tax.speciesinfo[tax.typelabel]
end

function _calcsimilarity(::Taxonomy, ::AbstractArray)
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
        M <: AbstractMatrix}(zmatrix::M, names::Vector{String})
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

function _getnames(gt::GeneralTypes, ::Bool)
    return gt.names
end

function _calcsimilarity(gt::GeneralTypes, ::AbstractArray)
    return gt.z
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

Metacommunity(abundances::AbstractArray, part::AbstractPartition, types::AbstractTypes)

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
    sourceabundances::A
    internalabundances::A
    types::Sim
    partition::Part
    ordinariness::Nullable{A}

    function (::Type{Metacommunity{FP, A, Sim, Part}}){FP <: AbstractFloat,
        A <: AbstractArray,
        Sim <: AbstractTypes,
        Part <: AbstractPartition}(abundances::A, types::Sim, part::Part)
        mcmatch(abundances, types, part) ||
            throw(ErrorException("Type or size mismatch between abundance array, partition and type list"))
        internalabundances = calcabundance(types, abundances)
        new{FP, A, Sim, Part}(abundances, internalabundances, types, part, Nullable{A}())
    end

    function (::Type{Metacommunity{FP, A, Sim, Part}}){FP <: AbstractFloat,
        A <: AbstractArray,
        Sim <: AbstractTypes,
        Part <: AbstractPartition}(abundances::A,
                                   meta::Metacommunity{FP, A, Sim, Part})
        mcmatch(abundances, meta.types, meta.part) ||
            throw(ErrorException("Type or size mismatch between abundance array, partition and type list"))
        internalabundances = calcabundance(types, abundances)
        new{FP, A, Sim, Part}(abundances, internalabundances, meta.types,
                              meta.part, meta.ordinariness)
    end
end

function Metacommunity{A <: AbstractArray,
    Meta <: AbstractMetacommunity}(abundances::A, meta::Meta)
    types = gettypes(meta)
    part = getpartition(meta)
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(A), A,
        typeof(types), typeof(part)}(abundances, types, part)
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
        typeof(types), typeof(part)}(ab, types, part)
    end
end

function Metacommunity{M <: AbstractMatrix, Sim <: AbstractTypes,
    Part <: AbstractPartition}(abundances::M,
                               types::Sim = UniqueTypes(size(abundances, 1)),
                               part::Part = Subcommunities(size(abundances, 2)))
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(M), M,
        typeof(types), typeof(part)}(abundances, types, part)
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
        typeof(types), typeof(part)}(ab, types, part)
    end
end

function Metacommunity{V <: AbstractVector, Sim <: AbstractTypes,
    Part <: AbstractPartition}(abundances::V,
                               types::Sim = UniqueTypes(size(abundances, 1)),
                               part::Part = Onecommunity())
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(V), V,
        typeof(types), typeof(part)}(abundances, types, part)
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
        typeof(types), typeof(part)}(ab, types, part)
    end
end

function Metacommunity{V <: AbstractVector,
    M <: AbstractMatrix}(abundances::V, zmatrix::M)
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(V), V,
        GeneralTypes, Onecommunity}(abundances,
                                    GeneralTypes(zmatrix),
                                    Onecommunity())
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
        GeneralTypes, Onecommunity}(ab, GeneralTypes(zmatrix), Onecommunity())
    end
    
end

function Metacommunity{M <: AbstractMatrix}(abundances::M, zmatrix::M)
    if sum(abundances) ≈ one(eltype(abundances))
        return Metacommunity{eltype(M), M,
        GeneralTypes, Subcommunities}(abundances,
                                      GeneralTypes(zmatrix),
                                      Subcommunities(size(abundances, 2)))
    else
        warn("Abundances not normalised to 1, correcting...")
        ab = abundances / sum(abundances)
        return Metacommunity{eltype(ab), typeof(ab),
        GeneralTypes, Subcommunities}(ab, GeneralTypes(zmatrix),
                                      Subcommunities(size(ab, 2)))
    end
end

function _gettypes{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.types
end

function _getpartition{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part})
    return meta.partition
end

function _getabundance{FP, A, Sim, Part}(meta::Metacommunity{FP, A, Sim, Part},
                                         input::Bool)
    return input ? meta.sourceabundances : meta.internalabundances
end

function _getordinariness!(meta::Metacommunity)
    if isnull(meta.ordinariness)
        meta.ordinariness = _calcordinariness(meta.types, meta.sourceabundances)
    end
    get(meta.ordinariness)
end
