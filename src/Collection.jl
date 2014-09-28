@doc """
### Abstract Similarity supertype

This type is the abstract superclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.
""" ->
abstract Similarity

@doc """
### Unique type, a subtype of Similarity

This type is the simplest Similarity subtype, which identifies all
individuals as unique and completely distinct from each other.
""" ->
immutable Unique <: Similarity
end

typealias Species Unique

@doc """
### Taxonomy type, a subtype of Similarity

This subtype of Similarity allows taxonomic similarity matrices
""" ->
immutable Taxonomy <: Similarity
    labels::Dict{String, (Float64, Dict{String, String})}
end

@doc """
### GeneralSimilarity, a general matrix-based Similarity subtype

This subtype of Similarity simply holds a matrix with similarities
between individuals.

### Members

**matrix** A two-dimensional matrix representing similarity between
           individuals. By default this will be the identity matrix,
           but will require the number of species to be instantiated.
""" ->
immutable GeneralSimilarity <: Similarity
    matrix::Matrix

    function GeneralSimilarity(z::Array{Float64, 2})
        size(z, 1) == size(z, 2) ||
        error("Similarity matrix is not square")
        
        isapprox(max(maximum(z), 1.), 1.) ||
        warn("Similarity matrix has values above 1")

        isapprox(min(minimum(z), 0.), 0.) ||
        warn("Similarity matrix has values below 0")

        new(z)
    end
end

getsimilarities{FP}(abundances::Array{FP}, s::Unique) =
    convert(Array{FP, 2}, eye(size(abundances, 1)))
    
getsimilarities{FP}(abundances::Array{FP}, s::Taxonomy) =
    error("Can't generate a taxonomic similarity matrix yet")

getsimilarities{FP}(abundances::Array{FP}, s::GeneralSimilarity) =
    convert(Array{FP, 2}, matrix)

@doc """
### Abstract Partition supertype

This type is the abstract superclass of all partitioning types.
Partition subtypes allow you to define how to partition your total
collection (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
""" ->
abstract Partition

immutable Subcommunity <: Partition; end
immutable Onecommunity <: Partition; end

legalpartition(abundances::Array, p::Subcommunity) =
    (ndims(abundances) == 2)

legalpartition(abundances::Array, p::Onecommunity) =
    (ndims(abundances) == 1)

@doc """
### Collection type

Type representing a single community or collection of communities. It
contains a collection of individuals which *may* be further
partitioned into smaller groups. For instance this may be an
ecosystem, which consists of a series of subcommunities.

The type stores relative abundances of different types, e.g. species,
and also allows for similarity between individuals.

### Parameterisation

**Collection{S, P, FP}**

**S** is the similarity type, e.g. Species, a subtype of Similarity.

**P** is the partition type, e.g. Subcommunity, a subtype of Partition.

**FP** is the kind of number storage, a subtype of FloatingPoint.

### Members

**abundances** An array of relative abundances. The first dimension
               represents the species, and further dimensions
               represent the structure of collection.

**Z** A two-dimensional matrix representing similarity between
      individuals of the base type, S. By default this will be the
      identity matrix.
""" ->
type Collection{S <: Similarity,
                P <: Partition,
                FP <: FloatingPoint}
    similarity::S
    partition::P
    abundances::Array{FP}
    similarities::Matrix{FP}
    normalise::Bool

    function Collection(s, p, ab, bool)
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = bool ? (ab / sum(ab)) : ab
        isapprox(sum(abundances), 1.0) || warn("Not normalised")
        sim = getsimilarities(abundances, partition)
        new(s, p, abundances, sim, bool)
    end

    function Collection(s, p, ab)
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = getsimilarities(abundances, partition)
        new(s, p, abundances, sim, true)
    end
    
    function Collection(s::S, ab)
        p = P()
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = getsimilarities(abundances, partition)
        new(s, p, abundances, sim, true)
    end
    
    function Collection(p::P, ab)
        legalpartition(ab, p) || error("Not a legal partition")        
        abundances = (ab / sum(ab))
        sim = getsimilarities(abundances, partition)
        new(S(), p, abundances, sim, true)
    end
    
    function Collection(ab)
        p = P()        
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = (ab / sum(ab))
        sim = getsimilarities(abundances, partition)
        new(S(), p, abundances, sim, true)
    end
end

Ecosystem{S, FP}(s::S, ab::Array{FP}) = Collection{S, Subcommunity, FP}(s, ab)
Ecosystem{FP}(ab::Array{FP}) = Collection{Species, Subcommunity, FP}(ab)

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



