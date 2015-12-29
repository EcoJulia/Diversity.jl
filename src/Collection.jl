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

- `matrix` A two-dimensional matrix representing similarity between
           individuals. By default this will be the identity matrix,
           but will require the number of species to be instantiated.
"""
type GeneralSimilarity{S <: AbstractFloat} <: Similarity
    
    """
    A two-dimensional matrix representing similarity between
    individuals. By default this will be the identity matrix,
    but will require the number of species to be instantiated.
    """
    z::Matrix{S}
    
    """
    ### Constructor for GeneralSimilarity
    
    Creates an instance of the GeneralSimilarity class, with an arbitrary similarity matrix.
    
    #### Arguments:
    - `z`: similarity matrix
    """
    function GeneralSimilarity(z::Matrix{S})
        size(z, 1) == size(z, 2) ||
        error("Similarity matrix is not square")
        
        isapprox(min(minimum(z), 0.), 0.) ||
        error("Similarity matrix has values below 0")
        
        isapprox(max(maximum(z), 1.), 1.) ||
        warn("Similarity matrix has values above 1")
        
        new(z)
    end
end

GeneralSimilarity{S <: AbstractFloat}(z::Matrix{S}) = GeneralSimilarity{S}(z)

get_similarities{S}(abundances::Array{S}, sim::Unique) =
    eye(Array(S, (size(abundances, 1), size(abundances, 1))))
    
get_similarities{S}(::Array{S}, sim::Taxonomy) =
    error("Can't generate a taxonomic similarity matrix yet")

get_similarities{S}(abundances::Array{S}, sim::GeneralSimilarity{S}) =
    size(abundances, 1) == size(sim.z, 1) ? sim.z : error("Similarity matrix size $(size(sim.z)) mismatch with number of types $(size(abundances, 1))")

"""
### Abstract Partition supertype for all partitioning types

This type is the abstract superclass of all partitioning types.
Partition subtypes allow you to define how to partition your total
collection (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
"""
abstract Partition

"""
### Partition type with multiple subccomunities
"""
immutable Subcommunity <: Partition; end

"""
### Partition type allowing only one subcommunity
"""
immutable Onecommunity <: Partition; end

legalpartition(abundances::Array, p::Subcommunity) =
    (ndims(abundances) == 2)

legalpartition(abundances::Array, p::Onecommunity) =
    (ndims(abundances) == 1)

"""
### Collection type, representing a collection of one or more subcommunities

Type representing a single community or collection of communities. It
contains a collection of individuals which *may* be further
partitioned into smaller groups. For instance this may be an
ecosystem, which consists of a series of subcommunities.

The type stores relative abundances of different types, e.g. species,
and also allows for similarity between individuals.

#### Parameterisation:

**Collection{S, P, FP}**

- `S` is the similarity type, e.g. Species, a subtype of Similarity.

- `P` is the partition type, e.g. Subcommunity, a subtype of Partition.

- `FP` is the kind of number storage, a subtype of AbstractFloat.

#### Members:

- `abundances` An array of relative abundances. The first dimension
               represents the species, and further dimensions
               represent the structure of collection.

- `z` A two-dimensional matrix representing similarity between
      individuals of the base type, S. By default this will be the
      identity matrix.
"""
type Collection{S <: Similarity,
                P <: Partition,
                FP <: AbstractFloat}
    similarity::S
    partition::P
    abundances::Array{FP}
    similarities::Matrix{FP}
    normalise::Bool

    function Collection(s, p, ab, normalise)
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = normalise ? (ab / sum(ab)) : ab
        isapprox(sum(abundances), 1.0) || warn("Not normalised")
        sim = get_similarities(abundances, partition)
        new(s, p, abundances, sim, normalise)
    end

    function Collection(s, p, ab)
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = get_similarities(abundances, partition)
        new(s, p, abundances, sim, true)
    end
    
    function Collection(s::S, ab)
        p = P()
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = get_similarities(abundances, partition)
        new(s, p, abundances, sim, true)
    end
    
    function Collection(p::P, ab)
        legalpartition(ab, p) || error("Not a legal partition")        
        abundances = (ab / sum(ab))
        sim = get_similarities(abundances, partition)
        new(S(), p, abundances, sim, true)
    end
    
    function Collection(ab)
        p = P()        
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = get_similarities(abundances, partition)
        new(S(), p, abundances, sim, true)
    end
end

"""
### Ecosystem type, representing an ecosystem of multiple subcommunities
"""
:Ecosystem

Ecosystem{S, FP}(s::S, ab::Array{FP}) = Collection{S, Subcommunity, FP}(s, ab)
Ecosystem{FP}(ab::Array{FP}) = Collection{Species, Subcommunity, FP}(ab)

"""
### Community type, representing a single community
"""
:Community

Community{S, FP}(s::S, ab::Array{FP}) = Collection{S, Onecommunity, FP}(s, ab)
Community{FP}(ab::Array{FP}) = Collection{Species, Onecommunity, FP}(ab)

getpartitions(c::Collection) = [1:ndims(c.abundances)]

getindividual(c::Collection) = -1

function summarydimensions(c::Collection, level::Integer)
    n = ndims(c.abundances)
    if level == 0
        return []
    elseif level == -1
        return [2:n]
    elseif level > 0 && level <= n
        return [1:level]
    end
end
