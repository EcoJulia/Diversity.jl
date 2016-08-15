"""
### AbstractPartition supertype for all partitioning types

This type is the abstract superclass of all partitioning types.
AbstractPartition subtypes allow you to define how to partition your total
supercommunity (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
"""
abstract AbstractPartition

"""
### Partition type with multiple subcommunities
"""
type Subcommunities <: AbstractPartition
    abundances::Matrix
    FPType::Type
    function Subcommunities{FP <: AbstractFloat}(abundances::Matrix{FP}, normalise::Bool = false)
        relative = normalise ? abundances / sum(abundances) : abundances
        isapprox(sum(relative), 1.0) || error("Not normalised")
        new(relative, FP)
    end
end

"""
### Partition type allowing only one subcommunity
"""
type Onecommunity <: AbstractPartition
    abundances::Vector
    FPType::Type
    function Onecommunity{FP <: AbstractFloat}(abundances::Vector{FP}, normalise::Bool = false)
        relative = normalise ? abundances / sum(abundances) : abundances
        isapprox(sum(relative), 1.0) || error("Not normalised")
        new(relative, FP)
    end
end

"""
### AbstractSimilarity supertype for all similarity measures

This type is the abstract superclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.
"""
abstract AbstractSimilarity

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
immutable Taxonomy <: AbstractSimilarity
    labels::Dict{AbstractString, Tuple{Float64, Dict{AbstractString, AbstractString}}}
end    

"""
### A general matrix-based AbstractSimilarity subtype

This subtype of AbstractSimilarity simply holds a matrix with similarities
between individuals.

#### Members:

- `zMatrix` A two-dimensional matrix representing similarity between
    individuals.
"""
type MatrixSimilarity <: AbstractSimilarity
    
    """
    A two-dimensional matrix representing similarity between
    individuals.
    """
    z::Matrix
    FPType::Type
    """
    ### Constructor for MatrixSimilarity
    
    Creates an instance of the MatrixSimilarity class, with an arbitrary similarity matrix.
    
    #### Arguments:
    - `z`: similarity matrix
    """
    function MatrixSimilarity{FP <: AbstractFloat}(z::Matrix{FP})
        size(z, 1) == size(z, 2) ||
        throw(DimensionMismatch("Similarity matrix is not square"))
        
        isapprox(min(minimum(z), 0.), 0.) ||
        throw(DomainError())
        
        isapprox(max(maximum(z), 1.), 1.) ||
        warn("Similarity matrix has values above 1")
        
        new(z, FP)
    end
end

function match!(part::AbstractPartition, ::Unique)
    true
end

function match!(part::AbstractPartition, ::Taxonomy)
    error("Taxonomic similarity not yet implemented")
end

function match!(part::AbstractPartition, sim::MatrixSimilarity)
    size(part.abundances, 1) == size(sim.z, 1) ||
    throw(DimensionMismatch("Similarity matrix size $(size(sim.z)) mismatch with number of types $(size(part.abundances, 1))"))
    if part.FPType != sim.FPType
        warn("Fixing Similarity type $(sim.FPType) to Partition type $(part.FPType)")
        sim.FPType = part.FPType
        sim.z = convert(Array{part.FPType}, sim.z)
    end
    true
end

"""
### AbstractSupercommunity supertype for all supercommunity types

This type is the abstract superclass of all supercommunity types.
AbstractSupercommunity subtypes allow you to define how to partition
your total supercommunity (e.g. an ecosystem) into smaller components
(e.g. subcommunities), and how to assess similarity between
individuals within it.
"""
abstract AbstractSupercommunity

"""
### Supercommunity type, representing a collection of individuals

Type representing a whole supercommunity containing a single community
or a collection of subcommunities. The supercommunity of individuals
*may* be further partitioned into smaller groups. For instance this
may be an ecosystem, which consists of a series of subcommunities. The
AbstractPartition subtype within it stores relative abundances of different
types, e.g. species, and also allows for similarity between
individuals.

#### Constructor:

**Supercommunity(part::AbstractPartition, sim::AbstractSimilarity)**


- `part` is an instance of type Part, the partition type, e.g.
  Subcommunities, a subtype of AbstractPartition.

- `sim` is an instance of type Sim, the similarity type, e.g. Species,
  a subtype of AbstractSimilarity.

#### Members:

- `partition` the instance of the AbstractPartition subtype, containing the
  subcommunities. These should be accessed through
  getabundance(::Supercommunity).

- `similarity` The instance of the AbstractSimilarity subtype, from which
  similarities between individuals can be calculated.

- `ordinariness` A cache of the ordinariness of the individuals in the
  Partition. Should only be accessed through
  getordinariness!(::Supercommunity), which will populate the cache if
  it has not yet been calculated.

- `FPType` is the kind of number storage, a subtype of AbstractFloat.

"""
type Supercommunity <: AbstractSupercommunity
    
    partition::AbstractPartition
    similarity::AbstractSimilarity
    ordinariness::Nullable{Array}
    FPType::Type

    function Supercommunity(part::AbstractPartition,
                            sim::AbstractSimilarity = Unique())
        match!(part, sim)
        new(part, sim, Nullable(), part.FPType)
    end
end

"""
### Ecosystem constructor for Supercommunity, representing an ecosystem of multiple subcommunities
"""
Ecosystem{FP <: AbstractFloat}(ab::Matrix{FP},
                               sim::AbstractSimilarity = Unique()) =
              Supercommunity(Subcommunities(ab), sim)

"""
### SingleCommunity contructor for Supercommunity, representing a single community
"""
SingleCommunity{FP <: AbstractFloat}(ab::Vector{FP},
                                     sim::AbstractSimilarity = Unique()) =
                    Supercommunity(Onecommunity(ab), sim)

"""
### Retrieves (and possibly calculates) the similarity matrix for a supercommunity
"""
function similaritymatrix
end

getsimilarity{Part <: AbstractPartition}(part::Part, ::Unique) =
    convert(Array{eltype(part.abundances), 2}, eye(size(part.abundances, 1)))
    
getsimilarity(part::AbstractPartition, ::Taxonomy) =
    error("Can't generate a taxonomic similarity matrix yet")

function getsimilarity(part::AbstractPartition, sim::MatrixSimilarity)
    match!(part, sim)
    return sim.z
end

getsimilarity{Sup <: AbstractSupercommunity}(sup::Sup) =
    getsimilarity(sup.partition, sup.similarity)

"""
### Retrieves (and possibly calculates) the relative abundances of a supercommunity
"""
function getabundance
end

getabundance{Sup <: AbstractSupercommunity}(sup::Sup) = sup.partition.abundances

"""
### Calculates the ordinarinesses of the subcommunities in a supercommunity
"""
function getordinariness
end

function getordinariness{Part <: AbstractPartition}(part::Part, ::Unique)
    return part.abundances
end
    
function getordinariness{Part <: AbstractPartition}(part::Part, ::Taxonomy)
    error("Can't generate a taxonomic similarity matrix yet")
end

function getOrdinariness!(part::AbstractPartition, sim::MatrixSimilarity)
    match!(part, sim)
    sim.z * part.abundances
end

"""
### Retrieves (and possibly calculates) the ordinarinesses of the subcommunities in a supercommunity
"""
function getordinariness!{Sup <: AbstractSupercommunity}(sup::Sup)
    if isnull(sup.ordinariness)
        sup.ordinariness = getordinariness(sup.partition, sup.similarity)
    end
    get(sup.ordinariness)
end

"""
### Retrieves (and possibly calculates) the ordinarinesses of a whole supercommunity
"""
function getsuperordinariness!{Sup <: AbstractSupercommunity}(sup::Sup)
    ord = getordinariness!(sup)
    sumoversubcommunities(sup.partition, ord)
end

"""
### Retrieves (and possibly calculates) the relative weights of the subcommunities
"""
function getweight{Sup <: AbstractSupercommunity}(sup::Sup)
    ab = getabundance(sup)
    sumovertypes(sup.partition, ab)
end

## Now create the functions for the iterator interface for Partitions
## and Supercommunities
import Base.start, Base.next, Base.done, Base.eltype, Base.length

function start(::Onecommunity)
    (1,)
end

function next(::Onecommunity, state::Tuple{Int64})
    index, = state
    (one.abundances, (index + 1, ))
end

function done(one::Onecommunity, state::Tuple{Int64})
    index, = state
    index != 1
end

function eltype(one::Onecommunity)
    Vector{one.FPType}
end

function length(::Onecommunity)
    1
end

function start(::Subcommunities)
    (1,)
end

function next(sub::Subcommunities, state::Tuple{Int64})
    index, = state
    (sub.abundances[:, index], (index + 1,))
end

function done(sub::Subcommunities, state::Tuple{Int64})
    index, = state
    index > length(sub)
end

function eltype(sub::Subcommunities)
    Vector{sub.FPType}
end

function length(sub::Subcommunities)
    size(sub.abundances, 2)
end

function start(sup::AbstractSupercommunity)
    (1, start(sup.partition))
end

function next(sup::AbstractSupercommunity, state::Tuple{Int64, Tuple})
    index_sup, index_part = state
    item_part, index_part = next(sup.partition, index_part)
    item_sup = getordinariness!(sup)[:, index_sup]
    ((item_sup, item_part), (index_sup + 1, index_part))
end

function done(sup::AbstractSupercommunity, state::Tuple{Int64, Tuple})
    index_sup, state_part = state
    done(sup.partition, state_part)
end

function eltype(sup::AbstractSupercommunity)
    (Vector{sup.FPType}, eltype(sup.partition))
end

function length(sup::AbstractSupercommunity)
    length(sup.partition)
end

