"""
!!summary(Abstract Similarity supertype for all similarity measures)

This type is the abstract superclass of all similarity types. Its
subtypes allow you to define how similarity is measured between
individuals.
"""
abstract Similarity

"""
!!summary(A subtype of Similarity where all individuals are completely distinct)

This type is the simplest Similarity subtype, which identifies all
individuals as unique and completely distinct from each other.
"""
immutable Unique <: Similarity
end

"""
!!summary(A subtype of Similarity where all species are completely distinct)

This type is the simplest Similarity subtype, which identifies all
species as unique and completely distinct from each other.
"""
typealias Species Unique

"""
!!summary(A subtype of Similarity with similarity between related taxa)

This subtype of Similarity allows taxonomic similarity matrices
"""
:Taxonomy

if VERSION < v"0.4-"
    immutable Taxonomy <: Similarity
        labels::Dict{String, (Float64, Dict{String, String})}
    end
else
    immutable Taxonomy <: Similarity
        labels::Dict{String, Tuple{Float64, Dict{String, String}}}
    end    
end


"""
!!summary(A general matrix-based Similarity subtype)

This subtype of Similarity simply holds a matrix with similarities
between individuals.

#### Members:

**matrix** A two-dimensional matrix representing similarity between
           individuals. By default this will be the identity matrix,
           but will require the number of species to be instantiated.
"""
immutable GeneralSimilarity <: Similarity
    matrix::Matrix

"""
!!summary(Constructor for GeneralSimilarity)

Creates an instance of the GeneralSimilarity class, with an arbitrary similarity matrix.

#### Arguments:
* Z: similarity matrix
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
    convert(Array{FP, 2}, matrix)

"""
!!summary(Abstract Partition supertype for all partitioning types)

This type is the abstract superclass of all partitioning types.
Partition subtypes allow you to define how to partition your total
collection (e.g. an ecosystem) into smaller components (e.g.
subcommunities).
"""
abstract Partition

"""
!!summary(Partition type with multiple subccomunities)
"""
immutable Subcommunity <: Partition; end

"""
!!summary(Partition type allowing only one subcommunity)
"""
immutable Onecommunity <: Partition; end

legalpartition(abundances::Array, p::Subcommunity) =
    (ndims(abundances) == 2)

legalpartition(abundances::Array, p::Onecommunity) =
    (ndims(abundances) == 1)

"""
!!summary(Collection type, representing a collection of one or more subcommunities)

Type representing a single community or collection of communities. It
contains a collection of individuals which *may* be further
partitioned into smaller groups. For instance this may be an
ecosystem, which consists of a series of subcommunities.

The type stores relative abundances of different types, e.g. species,
and also allows for similarity between individuals.

#### Parameterisation:

**Collection{S, P, FP}**

**S** is the similarity type, e.g. Species, a subtype of Similarity.

**P** is the partition type, e.g. Subcommunity, a subtype of Partition.

**FP** is the kind of number storage, a subtype of FloatingPoint.

#### Members:

**abundances** An array of relative abundances. The first dimension
               represents the species, and further dimensions
               represent the structure of collection.

**Z** A two-dimensional matrix representing similarity between
      individuals of the base type, S. By default this will be the
      identity matrix.
"""
type Collection{S <: Similarity,
                P <: Partition,
                FP <: FloatingPoint}
    similarity::S
    partition::P
    abundances::Array{FP}
    similarities::Matrix{FP}
    normalise::Bool

    function Collection(s, p, ab, normalise)
        legalpartition(ab, p) || error("Not a legal partition")
        abundances = normalise ? (ab / sum(ab)) : ab
        isapprox(sum(abundances), 1.0) || warn("Not normalised")
        sim = getsimilarities(abundances, partition)
        new(s, p, abundances, sim, normalise)
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

"""
!!summary(Ecosystem type, representing an ecosystem of multiple subcommunities)
"""
:Ecosystem

Ecosystem{S, FP}(s::S, ab::Array{FP}) = Collection{S, Subcommunity, FP}(s, ab)
Ecosystem{FP}(ab::Array{FP}) = Collection{Species, Subcommunity, FP}(ab)

"""
!!summary(Community type, representing a single community)
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



