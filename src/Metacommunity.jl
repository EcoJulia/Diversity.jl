"""
### AbstractPartition metatype for all partitioning types

This type is the abstract metaclass of all partitioning types.
AbstractPartition subtypes allow you to define how to partition your total
metacommunity (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
"""
abstract AbstractPartition{FP <: AbstractFloat, A <: AbstractArray}

"""
### Partition type with multiple subcommunities
"""
type Subcommunities{FP, M <: AbstractMatrix} <: AbstractPartition{FP, M}
    abundances::M

    Subcommunities(m::AbstractMatrix{FP}) = new(m)
end

function Subcommunities{FP <: AbstractFloat}(abundances::AbstractMatrix{FP},
                                             normalise::Bool = false)
    relative = normalise ? abundances / sum(abundances) : abundances
    sum(relative) ≈ 1.0 || error("Not normalised")
    Subcommunities{FP, typeof(relative)}(relative)
end

function Subcommunities{IT <: Integer}(abundances::AbstractMatrix{IT})
    relative = abundances / sum(abundances)
    Subcommunities{eltype(relative), typeof(relative)}(relative)
end

"""
### Partition type allowing only one subcommunity
"""
type Onecommunity{FP, V <: AbstractVector} <: AbstractPartition{FP, V}
    abundances::V

    Onecommunity(v::AbstractVector{FP}) = new(v)
end

function Onecommunity{FP <: AbstractFloat}(abundances::AbstractVector{FP},
                                           normalise::Bool = false)
    relative = normalise ? abundances / sum(abundances) : abundances
    sum(relative) ≈ 1.0 || error("Not normalised")
    Onecommunity{FP, typeof(relative)}(relative)
end

function Onecommunity{IT <: Integer}(abundances::AbstractVector{IT})
    relative = abundances / sum(abundances)
    Onecommunity{eltype(relative), typeof(relative)}(relative)
end

"""
### AbstractSimilarity metatype for all similarity measures

This type is the abstract metaclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.
"""
abstract AbstractSimilarity{FP <: AbstractFloat, A <: AbstractMatrix}

"""
### A subtype of AbstractSimilarity where all individuals are completely distinct

This type is the simplest AbstractSimilarity subtype, which identifies all
individuals as unique and completely distinct from each other.
"""
immutable Unique <: AbstractSimilarity
end

"""
### A subtype of AbstractSimilarity where all species are completely distinct

This type is the simplest AbstractSimilarity subtype, which identifies all
species as unique and completely distinct from each other.
"""
typealias Species Unique

"""
### A subtype of AbstractSimilarity with similarity between related taxa

This subtype of AbstractSimilarity allows taxonomic similarity matrices
"""
immutable Taxonomy{FP, STR <: AbstractString, D <: Dict} <: AbstractSimilarity{FP}
    labels::D
end

Taxonomy{FP <: AbstractFloat, STR <: AbstractString}(labels::Dict{STR, Tuple{FP, Dict{STR, STR}}}) = Taxonomy{FP, STR, typeof(labels)}(labels)

"""
### A general matrix-based AbstractSimilarity subtype

This subtype of AbstractSimilarity simply holds a matrix with similarities
between individuals.

#### Members:

- `zMatrix` A two-dimensional matrix representing similarity between
    individuals.
"""
type MatrixSimilarity{FP, M <: AbstractMatrix} <: AbstractSimilarity{FP, M}
    """
    A two-dimensional matrix representing similarity between
    individuals.
    """
    z::M
end

"""
### Constructor for MatrixSimilarity

Creates an instance of the MatrixSimilarity class, with an arbitrary similarity matrix.

#### Arguments:
- `z`: similarity matrix
"""
function MatrixSimilarity{FP <: AbstractFloat}(z::AbstractMatrix{FP})
    size(z, 1) == size(z, 2) ||
    throw(DimensionMismatch("Similarity matrix is not square"))

    min(minimum(z), 0.0) ≈ 0.0 || throw(DomainError())

    max(maximum(z), 1.0) ≈ 1.0 || warn("Similarity matrix has values above 1")

    MatrixSimilarity{FP, typeof(z)}(z)
end

function psmatch(part::AbstractPartition, ::Unique)
    true
end

function psmatch{FP}(part::AbstractPartition{FP}, ::Taxonomy{FP})
    error("Taxonomic similarity not yet implemented")
end

function psmatch{FP}(part::AbstractPartition{FP}, sim::MatrixSimilarity{FP})
    size(part.abundances, 1) == size(sim.z, 1) ||
    throw(DimensionMismatch("Similarity matrix size $(size(sim.z)) mismatch with number of types $(size(part.abundances, 1))"))
    eltype(part.abundances) == eltype(sim.z) ||
    throw("Similarity type $(eltype(part.abundances)) does not match Partition type $(eltype(sim.z))")
    true
end

"""
### AbstractMetacommunity metatype for all metacommunity types

This type is the abstract metaclass of all metacommunity types.
AbstractMetacommunity subtypes allow you to define how to partition
your total metacommunity (e.g. an ecosystem) into smaller components
(e.g. subcommunities), and how to assess similarity between
individuals within it.
"""
abstract AbstractMetacommunity{FP <: AbstractFloat, A <: AbstractArray, Part <: AbstractPartition, Sim <: AbstractSimilarity}

"""
### Metacommunity type, representing a collection of individuals

Type representing a whole metacommunity containing a single community
or a collection of subcommunities. The metacommunity of individuals
*may* be further partitioned into smaller groups. For instance this
may be an ecosystem, which consists of a series of subcommunities. The
AbstractPartition subtype within it stores relative abundances of different
types, e.g. species, and also allows for similarity between
individuals.

#### Constructor:

**Metacommunity(part::AbstractPartition, sim::AbstractSimilarity)**


- `part` is an instance of type Part, the partition type, e.g.
  Subcommunities, a subtype of AbstractPartition.

- `sim` is an instance of type Sim, the similarity type, e.g. Species,
  a subtype of AbstractSimilarity.

#### Members:

- `partition` the instance of the AbstractPartition subtype, containing the
  subcommunities. These should be accessed through
  getabundance(::Metacommunity).

- `similarity` The instance of the AbstractSimilarity subtype, from which
  similarities between individuals can be calculated.

- `ordinariness` A cache of the ordinariness of the individuals in the
  Partition. Should only be accessed through
  getordinariness!(::Metacommunity), which will populate the cache if
  it has not yet been calculated.

- `FPType` is the kind of number storage, a subtype of AbstractFloat.

"""
type Metacommunity{FP, A, Part, Sim} <: AbstractMetacommunity{FP, A, Part, Sim}

    partition::Part
    similarity::Sim
    ordinariness::Nullable{A}
end

function Metacommunity{Part <: AbstractPartition, Sim <: AbstractSimilarity}(part::Part, sim::Sim)
    psmatch(part, sim) || throw("Type mismatch")
    A = typeof(part.abundances)
    FP = eltype(A)
    Metacommunity{FP, A, Part, Sim}(part, sim, Nullable{A}())
end

function Metacommunity{Part <: AbstractPartition}(part::Part)
    A = typeof(part.abundances)
    FP = eltype(A)
    Metacommunity{FP, A, Part, Unique}(part, Unique(), Nullable{A}())
end

"""
### Multiple subcommunity constructors for Metacommunity
"""
Metacommunity{Mat <: AbstractMatrix,
Sim <: AbstractSimilarity}(ab::Mat, sim::Sim = Unique()) =
    Metacommunity(Subcommunities(ab), sim)
Metacommunity{Mat1 <: AbstractMatrix, Mat2 <: AbstractMatrix}(ab::Mat1,
                                                               z::Mat2) =
    Metacommunity(Subcommunities(ab), MatrixSimilarity(z))

"""
### Single subcommunity contructor for Metacommunity
"""
Metacommunity{Vec <: AbstractVector,
Sim <: AbstractSimilarity}(ab::Vec, sim::Sim = Unique()) =
    Metacommunity(Onecommunity(ab), sim)

"""
### Retrieves (and possibly calculates) the similarity matrix for a metacommunity
"""
function similaritymatrix
end

getsimilarity{Part <: AbstractPartition}(part::Part, ::Unique) =
    convert(Array{eltype(part.abundances), 2}, eye(size(part.abundances, 1)))

getsimilarity(part::AbstractPartition, ::Taxonomy) =
    error("Can't generate a taxonomic similarity matrix yet")

function getsimilarity(part::AbstractPartition, sim::MatrixSimilarity)
    psmatch(part, sim)
    return sim.z
end

getsimilarity{Sup <: AbstractMetacommunity}(sup::Sup) =
    getsimilarity(sup.partition, sup.similarity)

"""
### Retrieves (and possibly calculates) the relative abundances of a metacommunity
"""
function getabundance
end

getabundance{Sup <: AbstractMetacommunity}(sup::Sup) = sup.partition.abundances

"""
### Calculates the ordinarinesses of the subcommunities in a metacommunity
"""
function getordinariness
end

# Generic get ordinariness method
function getordinariness{Part <: AbstractPartition,
    Sim <: AbstractSimilarity}(part::Part, sim::Sim)
    psmatch(part, sim)
    return sim.z * part.abundances
end

# Simpler for distinct types
function getordinariness{Part <: AbstractPartition}(part::Part, ::Unique)
    return part.abundances
end

"""
### Retrieves (and possibly calculates) the ordinarinesses of the subcommunities in a metacommunity
"""
function getordinariness!{Sup <: AbstractMetacommunity}(sup::Sup)
    if isnull(sup.ordinariness)
        sup.ordinariness = getordinariness(sup.partition, sup.similarity)
    end
    get(sup.ordinariness)
end

"""
### Retrieves (and possibly calculates) the ordinarinesses of a whole metacommunity
"""
function getmetaordinariness!{Sup <: AbstractMetacommunity}(sup::Sup)
    ord = getordinariness!(sup)
    sumoversubcommunities(sup.partition, ord)
end

"""
### Retrieves (and possibly calculates) the relative weights of the subcommunities
"""
function getweight{Sup <: AbstractMetacommunity}(sup::Sup)
    ab = getabundance(sup)
    sumovertypes(sup.partition, ab)
end

@inline function sumoversubcommunities{FP <: AbstractFloat,
    M <: AbstractMatrix}(sc::Subcommunities{FP, M}, vals::M)
    vec(mapslices(sum, vals, 2))::Vector{FP}
end

@inline function sumoversubcommunities{FP <: AbstractFloat,
    V <: AbstractVector}(oc::Onecommunity{FP, V}, vals::V)
    vals
end

@inline function sumovertypes{FP <: AbstractFloat,
    M <: AbstractMatrix}(sc::Subcommunities{FP, M}, vals::M)
    mapslices(sum, vals, 1)::M
end

@inline function sumovertypes{FP <: AbstractFloat,
    V <: AbstractVector}(oc::Onecommunity{FP, V}, vals::V)
    sum(vals)
end


## Now create the functions for the iterator interface for Partitions
## and Metacommunities
import Base.start, Base.next, Base.done, Base.eltype, Base.length

function start{FP <: AbstractFloat, V <: AbstractVector}(::Onecommunity{FP, V})
    (1,)
end

function next{FP <: AbstractFloat, V <: AbstractVector}(one::Onecommunity{FP, V},
                                                        state::Tuple{Int64})
    index, = state
    (one.abundances, (index + 1, ))
end

function done{FP <: AbstractFloat, V <: AbstractVector}(one::Onecommunity{FP, V},
                                                        state::Tuple{Int64})
    index, = state
    index != 1
end

function eltype{FP <: AbstractFloat, V <: AbstractVector}(one::Onecommunity{FP, V})
    V
end

function length{FP <: AbstractFloat, V <: AbstractVector}(::Onecommunity{FP, V})
    1
end

function start{FP <: AbstractFloat, M <: AbstractMatrix}(::Subcommunities{FP, M})
    (1,)
end

function next{FP <: AbstractFloat, M <: AbstractMatrix}(sub::Subcommunities{FP, M},
                                                        state::Tuple{Int64})
    index, = state
    (sub.abundances[:, index], (index + 1,))
end

function done{FP <: AbstractFloat, M <: AbstractMatrix}(sub::Subcommunities{FP, M},
                                                        state::Tuple{Int64})
    index, = state
    index > length(sub)
end

function eltype{FP <: AbstractFloat, M <: AbstractMatrix}(sub::Subcommunities{FP, M})
    Vector{FP}
end

function length{FP <: AbstractFloat, M <: AbstractMatrix}(sub::Subcommunities{FP, M})
    size(sub.abundances, 2)
end

function start{Sup <: AbstractMetacommunity}(sup::Sup)
    (1, start(sup.partition))
end

function next{Sup <: AbstractMetacommunity}(sup::Sup, state::Tuple{Int64, Tuple})
    index_sup, index_part = state
    item_part, index_part = next(sup.partition, index_part)
    item_sup = getordinariness!(sup)[:, index_sup]
    ((item_sup, item_part), (index_sup + 1, index_part))
end

function done{Sup <: AbstractMetacommunity}(sup::Sup, state::Tuple{Int64, Tuple})
    index_sup, state_part = state
    done(sup.partition, state_part)
end

function eltype{Sup <: AbstractMetacommunity}(sup::Sup)
    (Vector{eltype(sup.ordinariness)}, eltype(sup.partition))
end

function length{Sup <: AbstractMetacommunity}(sup::Sup)
    length(sup.partition)
end
