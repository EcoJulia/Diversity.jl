"""
### Abstract Partition supertype for all partitioning types

This type is the abstract superclass of all partitioning types.
Partition subtypes allow you to define how to partition your total
supercommunity (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
"""
abstract Partition

"""
### Partition type with multiple subcommunities
"""
type Subcommunities <: Partition
    abundances::Matrix
    FPType::Type
    function Subcommunities{FP <: AbstractFloat}(abundances::Matrix{FP}, normalise::Bool = false)
        relative = normalise ? abundances / sum(abundances) : abundances
        isapprox(sum(relative), 1.0) || warn("Not normalised")
        new(relative, FP)
    end
end

"""
### Partition type allowing only one subcommunity
"""
type Onecommunity <: Partition
    abundances::Vector
    FPType::Type
    function Onecommunity{FP <: AbstractFloat}(abundances::Vector{FP}, normalise::Bool = false)
        relative = normalise ? abundances / sum(abundances) : abundances
        isapprox(sum(relative), 1.0) || warn("Not normalised")
        new(relative, FP)
    end
end

"""
### Abstract Similarity supertype for all similarity measures

This type is the abstract superclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.
"""
abstract Similarity

"""
### A subtype of Similarity where all individuals are completely distinct

This type is the simplest Similarity subtype, which identifies all
individuals as unique and completely distinct from each other.
"""
immutable Unique <: Similarity
end

"""
### A subtype of Similarity where all species are completely distinct

This type is the simplest Similarity subtype, which identifies all
species as unique and completely distinct from each other.
"""
typealias Species Unique

"""
### A subtype of Similarity with similarity between related taxa

This subtype of Similarity allows taxonomic similarity matrices
"""
immutable Taxonomy <: Similarity
    labels::Dict{AbstractString, @compat(Tuple{Float64, Dict{AbstractString, AbstractString}})}
end    

"""
### A general matrix-based Similarity subtype

This subtype of Similarity simply holds a matrix with similarities
between individuals.

#### Members:

- `zMatrix` A two-dimensional matrix representing similarity between
    individuals.
"""
type MatrixSimilarity <: Similarity
    
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

function match!{Part <: Partition}(part::Part, ::Unique)
    true
end

function match!{Part <: Partition}(part::Part, ::Taxonomy)
    error("Taxonomic similarity not yet implemented")
end

function match!{Part <: Partition}(part::Part, sim::MatrixSimilarity)
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
### Supercommunity type, representing a collection of individuals

Type representing a whole supercommunity containing a single community
or a collection of subcommunities. The supercommunity of individuals
*may* be further partitioned into smaller groups. For instance this
may be an ecosystem, which consists of a series of subcommunities. The
Partition subtype within it stores relative abundances of different
types, e.g. species, and also allows for similarity between
individuals.

#### Parameterisation:

**Supercommunity{Sim, Part}**

- `Part` is the partition type, e.g. Subcommunities, a subtype of Partition.

- `Sim` is the similarity type, e.g. Species, a subtype of Similarity.

#### Members:

- `partition` the instance of the Partition subtype, containing the 

- `similarity` The instance of the Similarity subtype, from which
  similarities between individuals can be calculated.

- `ordinariness` A cache of the ordinariness of the individuals in the
  Partition. Should only be accessed through
  getOrdinariness!(::Supercommunity), which will populate the cache if
  it has not yet been calculated.

- `FPType` is the kind of number storage, a subtype of AbstractFloat.

"""
type Supercommunity
    
    partition::Partition
    similarity::Similarity
    ordinariness::Nullable{Array}
    FPType::Type

    function Supercommunity{Part <: Partition,
                            Sim <: Similarity}(part::Part, sim::Sim = Unique())
        match!(part, sim)
        new(part, sim, Nullable(), part.FPType)
    end
end

"""
### Ecosystem type, representing an ecosystem of multiple subcommunities
"""
:Ecosystem

Ecosystem{FP <: AbstractFloat,
          Sim <: Similarity}(ab::Matrix{FP}, sim::Sim = Unique()) =
              Supercommunity(Subcommunities(ab), sim)

"""
### Community type, representing a single community
"""
:Community

Community{FP <: AbstractFloat,
          Sim <: Similarity}(ab::Vector{FP}, sim::Sim = Unique()) =
              Supercommunity(Onecommunity(ab), sim)

getSimilarityMatrix{Part <: Partition}(part::Part, ::Unique) =
    convert(Array{part.FPType}, eye(size(part.abundances, 1)))
    
getSimilarityMatrix{Part <: Partition}(part::Part, ::Taxonomy) =
    error("Can't generate a taxonomic similarity matrix yet")

getSimilarityMatrix{Part <: Partition}(part::Part, sim::MatrixSimilarity) =
    match!(part, sim) && sim.z

getSimilarityMatrix{Sup <: Supercommunity}(sup::Sup) =
    getSimilarityMatrix(sup.partition, sup.similarity)

getAbundances{Sup <: Supercommunity}(sup::Sup) = sup.partition.abundances

function getOrdinariness!{Part <: Partition}(part::Part, ::Unique)
    part.abundances
end
    
function getOrdinariness!{Part <: Partition}(part::Part, ::Taxonomy)
    error("Can't generate a taxonomic similarity matrix yet")
end

function getOrdinariness!{Part <: Partition}(part::Part, sim::MatrixSimilarity)
    match!(part, sim)
    sim.z * part.abundances
end

function getOrdinariness!{Sup <: Supercommunity}(sup::Sup)
    isnull(sup.ordinariness) &&
    (sup.ordinariness = getOrdinariness!(sup.partition, sup.similarity))
    get(sup.ordinariness)
end




