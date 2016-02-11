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

if VERSION < v"0.4-"
@doc """
### A subtype of Similarity with similarity between related taxa

This subtype of Similarity allows taxonomic similarity matrices
""" ->
    immutable Taxonomy <: Similarity
        labels::Dict{String, (Float64, Dict{String, String})}
    end
else
@doc """
### A subtype of Similarity with similarity between related taxa

This subtype of Similarity allows taxonomic similarity matrices
""" ->
    immutable Taxonomy <: Similarity
        labels::Dict{AbstractString, Tuple{Float64, Dict{AbstractString, AbstractString}}}
    end    
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
immutable GeneralSimilarity <: Similarity
"""A two-dimensional matrix representing similarity between
individuals. By default this will be the identity matrix,
but will require the number of species to be instantiated."""
    matrix::Matrix

"""
### Constructor for GeneralSimilarity

Creates an instance of the GeneralSimilarity class, with an arbitrary similarity matrix.

#### Arguments:
- `Z`: similarity matrix
"""
    function GeneralSimilarity(Z::Array{Float64, 2})
        size(Z, 1) == size(Z, 2) ||
        error("Similarity matrix is not square")
        
        isapprox(max(maximum(Z), 1.), 1.) ||
        warn("Similarity matrix has values above 1")

        isapprox(min(minimum(Z), 0.), 0.) ||
        warn("Similarity matrix has values below 0")

        new(Z)
    end
end

getsimilarities{FP}(abundances::Array{FP}, s::Unique) =
    convert(Array{FP, 2}, eye(size(abundances, 1)))
    
getsimilarities{FP}(abundances::Array{FP}, s::Taxonomy) =
    error("Can't generate a taxonomic similarity matrix yet")

getsimilarities{FP}(abundances::Array{FP}, s::GeneralSimilarity) =
    convert(Array{FP, 2}, s.matrix)

"""
### Abstract Supercommunity supertype for all partitioning types

This type is the abstract superclass of all partitioning types.
Supercommunity subtypes allow you to define how to partition your total
supercommunity (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
"""
abstract Supercommunity

"""
### Supercommunity type with multiple subccomunities
"""
immutable Subcommunities <: Supercommunity; end

"""
### Supercommunity type allowing only one subcommunity
"""
immutable Onecommunity <: Supercommunity; end

legalpartition(abundances::Array, p::Subcommunities) =
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

**Collection{S, Super, FP}**

- `Sim` is the similarity type, e.g. Species, a subtype of Similarity.

- `Super` is the supercommunity type, e.g. Subcommunities, a subtype of Supercommunity.

- `FP` is the kind of number storage, a subtype of AbstractFloat.

#### Members:

- `abundances` An array of relative abundances. The first dimension
               represents the species, and further dimensions
               represent the structure of collection.

- `Z` A two-dimensional matrix representing similarity between
      individuals of the base type, Sim. By default this will be the
      identity matrix.
"""
type Collection{Sim <: Similarity,
                Super <: Supercommunity,
                FP <: AbstractFloat}
    similarity::Sim
    supercommunity::Super
    abundances::Array{FP}
    similarities::Matrix{FP}
    normalise::Bool

    function Collection(s, sup, ab, normalise)
        legalpartition(ab, sup) || error("Not a legal partition")
        abundances = normalise ? (ab / sum(ab)) : ab
        isapprox(sum(abundances), 1.0) || warn("Not normalised")
        sim = getsimilarities(abundances, s)
        new(s, sup, abundances, sim, normalise)
    end

    function Collection(s, sup, ab)
        legalpartition(ab, sup) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = getsimilarities(abundances, s)
        new(s, sup, abundances, sim, true)
    end
    
    function Collection(s::Sim, ab)
        sup = Super()
        legalpartition(ab, sup) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = getsimilarities(abundances, s)
        new(s, sup, abundances, sim, true)
    end
    
    function Collection(sup::Super, ab)
        legalpartition(ab, sup) || error("Not a legal partition")        
        abundances = (ab / sum(ab))
        s = Sim()
        sim = getsimilarities(abundances, s)
        new(s, sup, abundances, sim, true)
    end
    
    function Collection(ab)
        sup = Super()        
        legalpartition(ab, sup) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        s = Sim()
        sim = getsimilarities(abundances, s)
        new(s, sup, abundances, sim, true)
    end
end

"""
### Ecosystem type, representing an ecosystem of multiple subcommunities
"""
typealias Ecosystem{Sim, FP} Collection{Sim, Subcommunities, FP}
